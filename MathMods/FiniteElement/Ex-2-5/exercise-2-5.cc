/* $Id: step-3.cc 24232 2011-09-02 09:47:37Z kronbichler $ */
/* Author: Wolfgang Bangerth, 1999, Guido Kanschat, 2011 */

/*    $Id: step-3.cc 24232 2011-09-02 09:47:37Z kronbichler $       */
/*                                                                */
/*    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2006, 2007, 2008, 2010, 2011 by the deal.II authors */
/*                                                                */
/*    This file is subject to QPL and may not be  distributed     */
/*    without copyright and license information. Please refer     */
/*    to the file deal.II/doc/license.html for the  text  and     */
/*    further information on this license.                        */

/* Modified for Exercise 2.5 of the finite element lecture 
   in Hamburg in Summer 2014 by W. Wollner
*/

#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/numerics/data_out.h>
#include <fstream>
#include <iostream>

using namespace dealii;


class Solution : public Function<2>
{
public:
  Solution () : Function<2>() {}
  
  double value (const Point<2>   &p,
		const unsigned int  component = 0) const;
  
  Tensor<1,2> gradient (const Point<2>   &p,
			const unsigned int  component = 0) const;
};

double Solution::value(const Point<2>   &p, const unsigned int) const
{
  return sin(M_PI * p(0))*sin(2.*M_PI * p(1));
}

Tensor<1,2> Solution::gradient (const Point<2>   &p,
				const unsigned int) const
{
  //EXERCISE: This is used to evaluate the gradient of the 
  //known reference solution and needs to be implemented. 
  //If you don't know how to do this have a look into the 
  //deal.II steps 1-3 and the function evaluation 
  //the value of this function above.
  Tensor<1,2> return_value;
  return_value[0] = M_PI * cos(M_PI * p(0)) * sin(2.*M_PI * p(1));
  return_value[1] = 2 * M_PI * sin(M_PI * p(0)) * cos(2.*M_PI * p(1));
  return return_value;
}

class Problem
{
  public:
    Problem (unsigned int deg);
    void run ();
    void summarize_results () const;

  private:
    void refine_grid(); 
    void setup_system ();
    void assemble_system ();
    void solve ();
    void output_results ();
    virtual void make_grid (unsigned int ref);

    Triangulation<2>     triangulation;
   
    FE_Q<2>              fe;
    DoFHandler<2>        dof_handler;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Vector<double>       solution;
    Vector<double>       system_rhs;
    std::vector<double> dofs;
    std::vector<double> l2_value;
    std::vector<double> h1_value;
    unsigned int n_iter;
    unsigned int iter;
};

Problem::Problem (unsigned int deg)
		:
                fe (deg),
		dof_handler (triangulation)
{
  n_iter = 5;
  dofs.resize(n_iter);
  l2_value.resize(n_iter);
  h1_value.resize(n_iter);
  iter = 0;
}

void Problem::make_grid (unsigned int ref)
{
  GridGenerator::hyper_cube (triangulation);
  triangulation.refine_global (ref);
  std::cout << "Number of active cells: "
	    << triangulation.n_active_cells()
	    << std::endl;
  std::cout << "Total number of cells: "
	    << triangulation.n_cells()
	    << std::endl;
}
void Problem::refine_grid ()
{
  triangulation.refine_global (1);
  std::cout << std::endl;
  std::cout << "Refining the triangulation ..."<<std::endl;
  std::cout << "Number of active cells: "
	    << triangulation.n_active_cells()
	    << std::endl;
  std::cout << "Total number of cells: "
	    << triangulation.n_cells()
	    << std::endl;
}
void Problem::setup_system ()
{
  dof_handler.distribute_dofs (fe);
  std::cout << "Number of degrees of freedom: "
	    << dof_handler.n_dofs()
	    << std::endl;
  dofs[iter] = dof_handler.n_dofs();
  
  CompressedSparsityPattern c_sparsity(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, c_sparsity);
  sparsity_pattern.copy_from(c_sparsity);

  system_matrix.reinit (sparsity_pattern);
  solution.reinit (dof_handler.n_dofs());
  system_rhs.reinit (dof_handler.n_dofs());
}

void Problem::assemble_system ()
{
  QGauss<2>  quadrature_formula(2);
  FEValues<2> fe_values (fe, quadrature_formula,
			 update_quadrature_points | update_values | update_gradients | update_JxW_values);
  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.size();

  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       cell_rhs (dofs_per_cell);

  std::vector<unsigned int> local_dof_indices (dofs_per_cell);

  DoFHandler<2>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      fe_values.reinit (cell);
      cell_matrix = 0;
      cell_rhs = 0;
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	for (unsigned int j=0; j<dofs_per_cell; ++j)
	  for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
	    cell_matrix(i,j) += (fe_values.shape_grad (i, q_point) *
				 fe_values.shape_grad (j, q_point) *
				 fe_values.JxW (q_point));
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
	{

 	  //EXERCISE
 	  /* This should implement the required right hand side
 	   *  so you need to implement the integrand of 
 	   *  \int_\Omega f \phi \, dx 
 	   *  Since you use a quadrature formula, you only need to write 
 	   *  The value in the quadrature point (q_point)
 	   *  which you can access using 
 	   *  fe_values.quadrature_point(q_point)
 	   *  Currently the righthand side f = 1 is implemented.
 	   *  The value of the test function is accessible using 
 	   *  fe_values.shape_value (i, q_point)
 	   *
 	   *  Note: do not remove the term fe_values.JxW (q_point)
 	   *  it contains the quadrature weights!
 	   */
	  Point<2> c_point = fe_values.quadrature_point(q_point);
	  cell_rhs(i) += (fe_values.shape_value (i, q_point) *
			  5 * M_PI * M_PI * sin(M_PI * c_point(0)) * sin(2 * M_PI * c_point(1)) *
			  fe_values.JxW (q_point));
	}
      cell->get_dof_indices (local_dof_indices);

      for (unsigned int i=0; i<dofs_per_cell; ++i)
	for (unsigned int j=0; j<dofs_per_cell; ++j)
	  system_matrix.add (local_dof_indices[i],
			     local_dof_indices[j],
			     cell_matrix(i,j));

      for (unsigned int i=0; i<dofs_per_cell; ++i)
	system_rhs(local_dof_indices[i]) += cell_rhs(i);
    }

  std::map<unsigned int,double> boundary_values;
  VectorTools::interpolate_boundary_values (dof_handler,
					    0,
					    ZeroFunction<2>(),
					    boundary_values);
  MatrixTools::apply_boundary_values (boundary_values,
				      system_matrix,
				      solution,
				      system_rhs);
}

void Problem::solve ()
{
  SolverControl           solver_control (10000, 1e-12);
  SolverCG<>              solver (solver_control);

  solver.solve (system_matrix, solution, system_rhs,
		PreconditionIdentity());
}

void Problem::output_results () 
{
  DataOut<2> data_out;
  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (solution, "solution");
  data_out.build_patches ();
  std::ofstream output ("solution.gpl");
  data_out.write_gnuplot (output);  

  Vector<float> difference_per_cell (triangulation.n_active_cells());

  VectorTools::integrate_difference (dof_handler,
				     solution,
				     Solution(),
				     difference_per_cell,
				     QGauss<2>(3),
				     VectorTools::L2_norm);
  l2_value[iter] = difference_per_cell.l2_norm();


  std::cout << "L2 Error "
	    << l2_value[iter]
  << std::endl;
  //EXERCISE: Here you should take care to evaluate the 
  // H1 seminorm.
  VectorTools::integrate_difference (dof_handler,
				     solution,
				     Solution(),
				     difference_per_cell,
				     QGauss<2>(3),
				     VectorTools::H1_seminorm);
  h1_value[iter] = difference_per_cell.l2_norm();
  std::cout << "H1 Error: "
	    << h1_value[iter]
	    << std::endl;

}

void Problem::summarize_results () const
{
  std::cout<<"DOFS\tL2 Norm\t\tEOC\tH1 Norm\t\tEOC"<<std::endl;
  std::cout<<"--------------------------------------------------"<<std::endl;
  for(unsigned int i = 0; i < n_iter; i++)
  {
    double order_p = 0.;
    double order_m = 0.;
    if(i > 1)
    {
      //EXERCISE: 
      /* The following lines need to calculate the estimated order of 
       * convergence for the L^2- and H^1 Norm.
       * To do so we have stored the values of the L^2 error on the different 
       * triangulations h_1, h_2, h_3 in the array named " l2_value "
       * and those for the H1 seminorm in h1_value.
       */
      order_p = log(l2_value[i-1] / l2_value[i]) / log(2);
      order_m = log(h1_value[i-1] / h1_value[i]) / log(2);
      std::cout<<dofs[i]<<"\t"<<l2_value[i]<<"\t"<<order_p<<"\t"<<h1_value[i]<<"\t"<<order_m<<std::endl;
    }
    else
    {
      std::cout<<dofs[i]<<"\t"<<l2_value[i]<<"\t---\t"<<h1_value[i]<<"\t---"<<std::endl;
    }
    

  }
}

void Problem::run ()
{
  make_grid (3);
  for(;iter < n_iter; iter++)
  {
    setup_system();
    assemble_system ();
    solve ();
    output_results (); 
    refine_grid();
  }
}

int main ()
{
  Problem problem(1);
  problem.run ();
  Problem problem_2(2);
  problem_2.run ();

  std::cout<<"With Bilinear elements:"<<std::endl;
  problem.summarize_results();

  std::cout<<std::endl;
  std::cout<<"With Biquadratic elements::"<<std::endl;
  problem_2.summarize_results();

  return 0;
}
