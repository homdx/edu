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

/* Modified for Exercise 3.5 of the finite element lecture 
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
#include <deal.II/numerics/vectors.h>
#include <deal.II/numerics/matrices.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/solution_transfer.h>

#include <fstream>
#include <iostream>

using namespace dealii;


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
    Vector<double>       last_solution;
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
  n_iter = 7;
  dofs.resize(n_iter);
  l2_value.resize(n_iter);
  h1_value.resize(n_iter);
  iter = 0;
}

void Problem::make_grid (unsigned int ref)
{
  GridGenerator::hyper_L (triangulation,0,2);
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
  dealii::SolutionTransfer<2,Vector<double> > mesh_transfer(dof_handler);

  triangulation.set_all_refine_flags();
  triangulation.prepare_coarsening_and_refinement();
  mesh_transfer.prepare_for_pure_refinement();
  triangulation.execute_coarsening_and_refinement();

  std::cout << std::endl;
  std::cout << "Refining the triangulation ..."<<std::endl;
  std::cout << "Number of active cells: "
	    << triangulation.n_active_cells()
	    << std::endl;
  std::cout << "Total number of cells: "
	    << triangulation.n_cells()
	    << std::endl;

  dof_handler.distribute_dofs (fe);
  last_solution.reinit(dof_handler.n_dofs());

  mesh_transfer.refine_interpolate(solution, last_solution);

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
  QGauss<2>  quadrature_formula(3);
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
	  //EXERCISE: Change Rhs For the second part of the exercises
	  cell_rhs(i) += (fe_values.shape_value (i, q_point) *
			  5.*M_PI*M_PI*
			  sin(M_PI * fe_values.quadrature_point(q_point)(0)) *
			  sin(2. * M_PI * fe_values.quadrature_point(q_point)(1)) *
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
  std::ofstream output ("solution.vtk");
  data_out.write_vtk (output);  

  if(iter > 0)
  {
    last_solution -= solution; 
    Vector<float> difference_per_cell (triangulation.n_active_cells());

    VectorTools::integrate_difference (dof_handler,
				       last_solution,
				       ZeroFunction<2>(),
				       difference_per_cell,
				       QGauss<2>(4),
				       VectorTools::L2_norm);
    l2_value[iter] = difference_per_cell.l2_norm();
    
    
    std::cout << "L2 Error "
	      << l2_value[iter]
	      << std::endl;
    
    VectorTools::integrate_difference (dof_handler,
				       last_solution,
				       ZeroFunction<2>(),
				       difference_per_cell,
				       QGauss<2>(4),
				       VectorTools::H1_seminorm);
    h1_value[iter] = difference_per_cell.l2_norm();
    std::cout << "H1 Error: "
	      << h1_value[iter]
	      << std::endl;
    last_solution = 0.;
  }

}

void Problem::summarize_results () const
{
  std::cout<<"DOFS\tL2 Norm\t\tEOC\tH1 Norm\t\tEOC"<<std::endl;
  std::cout<<"--------------------------------------------------"<<std::endl;
  for(unsigned int i = 0; i < n_iter; i++)
  {
    double order_p = 0.;
    double order_m = 0.;
    if(i > 2)
    {
      order_p = log(fabs(l2_value[i-2]-l2_value[i-1])/fabs(l2_value[i-1]-l2_value[i]) )/log(2);
      order_m = log(fabs(h1_value[i-2]-h1_value[i-1])/fabs(h1_value[i-1]-h1_value[i]) )/log(2);
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
  make_grid (1);
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
  Problem problem2(2);
  problem2.run ();

 
  std::cout<<"With Bilinear elements:"<<std::endl;
  problem.summarize_results();
  std::cout<<"With Biquadratic elements:"<<std::endl;
  problem2.summarize_results();

  return 0;
}
