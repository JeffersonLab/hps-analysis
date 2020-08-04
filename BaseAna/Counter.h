//
//  Counter.h
//
//
#ifndef __Counter__  
#define __Counter__

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wconversion"
#include "TROOT.h"
#pragma clang diagnostic pop

#include <iostream>
#include <iomanip>

using namespace std;

class Counter : public TObject {
  
public:
  int pass;
  int fail;
  string name;
  
public:
  Counter(string n="test"){
    pass=0;
    fail=0;
    name = n;
  }
	Counter(Counter *c){
		pass= c->pass;
		fail= c->fail;
		name= c->name;
	}

  ~Counter(){};
  bool Test(bool t){ ///< Implements a Test(). Store a pass if true.
    if(t) pass++;
    else  fail++;
    return(t);
  }
	bool ATest(bool t){ ///< Implements an anti-Test, Store a pass if false.
    if(t) fail++;
    else  pass++;
    return(t);
  }

  void PrintOut(int num_evt) const{ ///< Print the number of passes and fails on one line, including a \% pass.
    
    cout << std::setw(25) << std::left << name << std::right
		<< " pass: " << std::setw(6) << pass
    << " = " << std::fixed << std::showpoint << std::setprecision(4) << std::setw(8)
		<< (100.*(double)pass)/(pass+fail)
    << "%  fail: " << std::setw(6) << fail
		<< " = " << std::setw(8) << (100.*(double)fail)/(pass+fail);
		if(num_evt>0) cout << "% Overall: "<< std::setw(8) << 100.*(double)pass/num_evt <<"%" << endl;
		else cout << endl;
  }
  
  void Print(Option_t *options="") const{ ///< Print the number of passes and fails on one line.
    TString opt = options;
    opt.ToLower();
    
    cout << std::setw(25) << std::left << name << std::right
		<< " pass: " << std::setw(6) << pass
    << " = " << std::fixed << std::showpoint << std::setprecision(4) << std::setw(8)
		<< (100.*(double)pass)/(pass+fail)
    << "%  fail: " << std::setw(6) << fail
		<< " = " << std::setw(8) << (100.*(double)fail)/(pass+fail);
		cout << endl;
  }
  
  ClassDef(Counter,1);
};



#endif
