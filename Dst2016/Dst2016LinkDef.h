#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ nestedclasses;

#pragma link C++ class Dst2016+;
#pragma link C++ class Multi_Value;
// #pragma link C++ class std::map<std::string, std::variant < vector<double>* , vector<int>* > >+;
#pragma link C++ class std::map<std::string, Multi_Value> >+;
#pragma link C++ class std::map<std::string, vector<double>* > +;
#pragma link C++ class std::map<std::string, vector<int>* > +;
#pragma link C++ class std::map<std::string, vector< vector<int> >* > +;
#endif
