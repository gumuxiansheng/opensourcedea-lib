
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
    
*/

package dea;

/**
 * A DEASolver exception. Used when the Lpsolve crashes. This should not happen has the check in the DEAProblem.solve() method
 * should throw an DEAException before bad data is sent to a model (e.g. throw MissingData DEAexception is there is no data in the
 * DEAProblem).
 * </br>
 * @author Hubert Virtos
 *
 */
public class DEASolverException extends DEAException {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	public DEASolverException() {
		super("The linear solver lpsolve encountered an error. This is likely caused by inconsistencies in the data sent to the Lpsolve.solveLPProblem method.");
	}
	
	public DEASolverException (String message) {
		super(message);
	}
}
