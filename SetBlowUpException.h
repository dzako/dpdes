/*  dPDEs - this program is an open research software performing rigorous integration in time of partial differential equations
    Copyright (C) 2010-2013  Jacek Cyranka

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
    
    Please consult the webpage www.cyranka.net,
    or contact me on jcyranka@gmail.com for further details.
*/

/*
 * SetBlowUpException.h
 *
 *  Created on: 22-01-2011
 *      Author: Jacek Cyranka
 */

#ifndef _CAPD_JACO_SETBLOWUPEXCEPTION_H_
#define _CAPD_JACO_SETBLOWUPEXCEPTION_H_

namespace capd{
namespace jaco{

///class used as a exception that is thrown when integrated set is too large, which result in blow up of the set
template<class ValueT>
class SetBlowUpException{
  typedef ValueT ValueType;
	ValueType m_value;
	const char* m_function;///name of a function that has thrown the exception
public:

	SetBlowUpException(const ValueType& value, const char* function) : m_value(value), m_function(function){}

	~SetBlowUpException(){
		//delete m_function;
	}

	friend std::ostream& operator<<(std::ostream& out, const SetBlowUpException& e) // output
	{
	  out<<"Set blow up occurred in ''"<<e.m_function<<"'' function, value that potentially caused blow up was "<<e.m_value<<
	      ". Probably one of the input parameters is not proper and caused blow up during the integration process,"<<
	      " for the usage refer the documentation, for the examples refer Section 7 in the paper."<<"\n";
	  return out;
	}


};

}}
#endif /* SETBLOWUPEXCEPTION_H_ */
