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
