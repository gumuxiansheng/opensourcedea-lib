
/*	<DEASolver (googleproject opensourcedea) is a java DEA solver.>
    Copyright (C) <2010>  <Hubert Paul Bernard VIRTOS>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
    
    
    @author Hubert Paul Bernard VIRTOS
    @version 0.1 2011-02-04
*/

package deaModels;

import dea.DEAProblem;
import dea.DEAVariableType;
import dea.SolverObjDirection;
import java.util.ArrayList;
import java.util.Arrays;
import linearSolver.*;



/**
 * The class implementing the CCR model.
 *
 *<p>
 *The linear problem of the Model is as follows:
 *<p>
 * <center>
 * <table border = "1">
 * <tr>
 * 		<td>Variable</td>
 * 		<td>Theta</td>
 * 		<td>Lambda 1</td>
 * 		<td>...</td>
 * 		<td>Lambda n</td>
 * 		<td>Slack 1</td>
 * 		<td>...</td>
 * 		<td>Slack p</td>
 * 		<td>DIR</td>
 * 		<td>RHS</td>
 * </tr>
 * <tr>
 * 		<td>Obj Coeff</td>
 * 		<td>1</td>
 * 		<td>0</td>
 * 		<td>0</td>
 * 		<td>0</td>
 * 		<td>0</td>
 * 		<td>0</td>
 * 		<td>0</td>
 * 		<td></td>
 * 		<td></td>
 * </tr>
 * <tr>
 * 		<td>Input 1</td>
 * 		<td>-Input 1,1</td>
 * 		<td>Input 1,1</td>
 * 		<td>...</td>
 * 		<td>Input 1, n</td>
 * 		<td>0</td>
 * 		<td>...</td>
 * 		<td>0</td>
 * 		<td>E</td>
 * 		<td>0</td>
 * </tr>
 * <tr>
 * 		<td>Input i</td>
 * 		<td>-Input i, 1</td>
 * 		<td>Input i, 1</td>
 * 		<td>...</td>
 * 		<td>Input i, n</td>
 * 		<td>0</td>
 * 		<td>...</td>
 * 		<td>0</td>
 * 		<td>E</td>
 * 		<td>0</td>
 * </tr>
 * <tr>
 * 		<td>Ouput p</td>
 * 		<td>0</td>
 * 		<td>Output p, 1</td>
 * 		<td>...</td>
 * 		<td>Output p, n</td>
 * 		<td>0</td>
 * 		<td>...</td>
 * 		<td>-1</td>
 * 		<td>E</td>
 * 		<td>Output p, 1</td>
 * </tr>
 * </table>
 * </center>
 * <p>
 *  Where the input values of the DMU being optimised are put in the Theta column and timed by -1
 *  (e.g. -Inputi 1 being the ith input of DMU1).
 * <p>
 * @author Hubert Virtos
 *
 */
public  class CCR {

	/**
	 * The method solving the CCR Problem
	 * @param deaP An instance of DEAProblem
	 */
	public static void solveCCR(DEAProblem deaP) {
		
		/* Declare & Collect the variables that will often be used in the process (rather
		 * than calling the different methods several times.*/
		int NbDMUs = deaP.getNumberOfDMUs();
		int NbVariables = deaP.getNumberOfVariables();
		double [] [] TransposedMatrix = new double [NbVariables] [NbDMUs];
		TransposedMatrix = deaP.getTranspose();
		
					
		
		
		/* As the CCR optimisation needs to be ran for all DMUs, 
		 * the program will loop through all DMUs.
		 * Also, the CCR model is solved in two phases.
		 * The problem will consequently be solved for each DMUs for Phase I and
		 * solved again for each DMUs for Phase II.*/
		
		  /////////////////////////////
		 //		Solve Phase I		//
		/////////////////////////////
		
		for (int i = 0; i < NbDMUs; i++) {
			
			//Build model for Phase I

			/*
			 * The model is built with an array as follows:
			 * 
			 * 	Variable:		Theta	  Lambda 1   .....	  Lambda n      Slack 1  ...  Slack p		DIR		RHS
			 *  
			 *  Obj Coeff:		  1		     0       .....       0		      0      ...     0			 
			 *  
			 *  DMU1 (input)   -Input1 1  Input1 1   .....     Input1 n       1      ...     0			 E		 0
			 *  DMUi (input)   -Inputi 1  Inputi 1   .....     Inputi n       0      ...     0			 E		 0
			 *  DMUp (output)     0	      Outputp 1  .....     Outputp n      0      ...    -1			 E	  Output p 1
			 * 
			 * 
			 *  Where the input values of the DMU being optimised are put in the Theta column and timed by -1
			 *  (e.g. -Inputi 1 being the ith input of DMU1).
			 */
			
			
			ArrayList<double[]> Constraints = new ArrayList<double []>();
			double[] ObjF = new double [NbDMUs + NbVariables + 1];
			double[] RHS1 = new double [NbVariables]; //RHS Phase I
			double[] RHS2 = new double [NbVariables]; //RHS Phase II

			for (int j = 0; j < NbVariables; j++) {
				
				//Build Model for each DMU
				
				//Build the Constraint Matrix
				double[] ConstraintRow = new double[NbDMUs + NbVariables + 1];
				if (deaP.getVariableType(j) == DEAVariableType.Input) {
					ConstraintRow[0] = TransposedMatrix[j] [i] * -1;
				}
				else {
					ConstraintRow[0] = 0;
				}
				System.arraycopy(TransposedMatrix[j], 0, ConstraintRow, 1, NbDMUs);
				if (deaP.getVariableType(j) == DEAVariableType.Input) {
					ConstraintRow[NbDMUs + 1 + j] = 1;
				}
				else {
					ConstraintRow[NbDMUs + 1 + j] = -1;
				}
				Constraints.add(ConstraintRow);
				
				
				//Build RHS
				if (deaP.getVariableType(j) == DEAVariableType.Input) {
					RHS1[j] = 0;
				}
				else {
					RHS1[j] = TransposedMatrix[j] [i];
				}
			
			}
		
			//Build Objective Function (Theta column is assigned the weight 1. All the other columns are left to 0).
			ObjF[0] = 1;
			
			
			//Solve
			SolverResults Sol = new SolverResults();
			Sol = Lpsolve.solveLPProblem(Constraints, ObjF, RHS1, SolverObjDirection.MIN);
			
			//Collect information from Phase I (Theta)
			deaP.setObjective(i, Sol.Objective);
			//deaP.Solution.Objectives[i] = Sol.Objective;
			deaP.setWeights(i, Sol.Weights);
			//deaP.Solution.Weights[i] = Sol.Weights;
			
			
			  /////////////////////////////
			 //		Solve Phase II		//
			/////////////////////////////
			/* The only things that are needed for Phase II are to:
			 * - add an extra Constraint Row to Constraints in order to ensure Theta Phase I
			 *   is not changed during Phase II Optimisation.
			 * - add the corresponding Theta to the RHS Array
			 * - change the Objective Function accordingly (all 1 on Slacks, all others coeff = 0).*/
			
			
			//Changing Constraint Matrix
			double[] ConstraintRow = new double[NbDMUs + NbVariables + 1];
			ConstraintRow[0] = 1;
			Constraints.add(ConstraintRow);
			
			//Changing RHS
			RHS2 = new double[NbVariables + 1];
			System.arraycopy(RHS1, 0, RHS2, 0, RHS1.length);
			RHS2[NbVariables] = Sol.Objective;
			
			//Change Objective Function
			Arrays.fill(ObjF,0);
			for (int j = NbDMUs + 1; j <= NbDMUs + NbVariables; j++) {
				ObjF[j] = 1;
			}
			
			//Solve the Phase II Problem
			Sol = Lpsolve.solveLPProblem(Constraints, ObjF, RHS2, SolverObjDirection.MAX);
			
			//Collect information from Phase II (Theta)
			System.arraycopy(Sol.VariableResult, 1, deaP.getLambdas(i)/*deaP._Solution.Lambdas[i]*/, 0, NbDMUs);
			System.arraycopy(Sol.VariableResult, NbDMUs + 1, deaP.getSlacks(i) /*deaP.Solution.Slacks[i]*/, 0, NbVariables);
			
			for (int j = 0; j < NbVariables; j++) {
				if(deaP.getVariableType(j) == DEAVariableType.Input) {
					//Projections
					deaP.setProjections(i, j, deaP.getObjective(i) * deaP.getDataMatrix(i, j) - deaP.getSlacks(i, j));
					/*deaP.Solution.Projections[i] [j] = 
						deaP.Solution.Objectives[i] * deaP.getDataMatrix(i, j) - deaP.Solution.Slacks[i] [j];*/
				}
				else {
					//Projections
					deaP.setProjections(i, j, deaP.getDataMatrix(i, j) + deaP.getSlacks(i, j));
					/*deaP.Solution.Projections[i] [j] =
						deaP.getDataMatrix(i, j) + deaP.Solution.Slacks[i] [j];*/
				}
			}
			
			
			
		}

		
	}
	
	
	

}
