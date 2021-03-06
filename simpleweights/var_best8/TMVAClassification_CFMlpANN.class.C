// Class: ReadCFMlpANN
// Automatically generated by MethodBase::MakeClass
//

/* configuration options =====================================================

#GEN -*-*-*-*-*-*-*-*-*-*-*- general info -*-*-*-*-*-*-*-*-*-*-*-

Method         : CFMlpANN::CFMlpANN
TMVA Release   : 4.0.7         [262151]
ROOT Release   : 5.27/06       [334598]
Creator        : wluo1
Date           : Tue Aug  7 22:38:25 2012
Host           : Linux lxbuild051.cern.ch 2.6.18-238.1.1.el5 #1 SMP Wed Jan 19 11:06:36 CET 2011 x86_64 x86_64 x86_64 GNU/Linux
Dir            : /afs/crc.nd.edu/user/w/wluo1/work/CMSSW_5_2_5/src/simpleMVA_ttbb
Training events: 1258
Analysis type  : [Classification]


#OPT -*-*-*-*-*-*-*-*-*-*-*-*- options -*-*-*-*-*-*-*-*-*-*-*-*-

# Set by User:
V: "False" [Verbose output (short form of "VerbosityLevel" below - overrides the latter one)]
H: "False" [Print method-specific help message]
NCycles: "2000" [Number of training cycles]
HiddenLayers: "N,N-1" [Specification of hidden layer architecture]
# Default:
VerbosityLevel: "Default" [Verbosity level]
VarTransform: "None" [List of variable transformations performed before training, e.g., "D_Background,P_Signal,G,N_AllClasses" for: "Decorrelation, PCA-transformation, Gaussianisation, Normalisation, each for the given class of events ('AllClasses' denotes all events of all classes, if no class indication is given, 'All' is assumed)"]
CreateMVAPdfs: "False" [Create PDFs for classifier outputs (signal and background)]
IgnoreNegWeightsInTraining: "False" [Events with negative weights are ignored in the training (but are included for testing and performance evaluation)]
##


#VAR -*-*-*-*-*-*-*-*-*-*-*-* variables *-*-*-*-*-*-*-*-*-*-*-*-

NVar 8
min_dr_tagged_jets            min_dr_tagged_jets            min_dr_tagged_jets            min_dr_tagged_jets                                              'F'    [0.504642367363,3.1229031086]
avg_dr_tagged_jets            avg_dr_tagged_jets            avg_dr_tagged_jets            avg_dr_tagged_jets                                              'F'    [0.824162423611,3.59096884727]
avg_tagged_dijet_mass         avg_tagged_dijet_mass         avg_tagged_dijet_mass         avg_tagged_dijet_mass                                           'F'    [45.5290412903,740.257873535]
second_dibjet_mass            second_dibjet_mass            second_dibjet_mass            second_dibjet_mass                                              'F'    [37.1900138855,661.453979492]
min_dr_jets                   min_dr_jets                   min_dr_jets                   min_dr_jets                                                     'F'    [0.504642367363,3.1229031086]
avg_dr_jets                   avg_dr_jets                   avg_dr_jets                   avg_dr_jets                                                     'F'    [0.870305836201,3.40389847755]
avg_dijet_mass                avg_dijet_mass                avg_dijet_mass                avg_dijet_mass                                                  'F'    [50.0762214661,740.257873535]
dijet_mass_second             dijet_mass_second             dijet_mass_second             dijet_mass_second                                               'F'    [49.8350563049,835.153930664]
NSpec 0


============================================================================ */

#include <vector>
#include <cmath>
#include <string>
#include <iostream>

#ifndef IClassifierReader__def
#define IClassifierReader__def

class IClassifierReader {

 public:

   // constructor
   IClassifierReader() : fStatusIsClean( true ) {}
   virtual ~IClassifierReader() {}

   // return classifier response
   virtual double GetMvaValue( const std::vector<double>& inputValues ) const = 0;

   // returns classifier status
   Bool_t IsStatusClean() const { return fStatusIsClean; }

 protected:

   Bool_t fStatusIsClean;
};

#endif

class ReadCFMlpANN : public IClassifierReader {

 public:

   // constructor
   ReadCFMlpANN( std::vector<std::string>& theInputVars ) 
      : IClassifierReader(),
        fClassName( "ReadCFMlpANN" ),
        fNvars( 8 ),
        fIsNormalised( true )
   {      
      // the training input variables
      const char* inputVars[] = { "min_dr_tagged_jets", "avg_dr_tagged_jets", "avg_tagged_dijet_mass", "second_dibjet_mass", "min_dr_jets", "avg_dr_jets", "avg_dijet_mass", "dijet_mass_second" };

      // sanity checks
      if (theInputVars.size() <= 0) {
         std::cout << "Problem in class \"" << fClassName << "\": empty input vector" << std::endl;
         fStatusIsClean = false;
      }

      if (theInputVars.size() != fNvars) {
         std::cout << "Problem in class \"" << fClassName << "\": mismatch in number of input values: "
                   << theInputVars.size() << " != " << fNvars << std::endl;
         fStatusIsClean = false;
      }

      // validate input variables
      for (size_t ivar = 0; ivar < theInputVars.size(); ivar++) {
         if (theInputVars[ivar] != inputVars[ivar]) {
            std::cout << "Problem in class \"" << fClassName << "\": mismatch in input variable names" << std::endl
                      << " for variable [" << ivar << "]: " << theInputVars[ivar].c_str() << " != " << inputVars[ivar] << std::endl;
            fStatusIsClean = false;
         }
      }

      // initialize min and max vectors (for normalisation)
      fVmin[0] = 0.504642367362976;
      fVmax[0] = 3.1229031085968;
      fVmin[1] = 0.824162423610687;
      fVmax[1] = 3.59096884727478;
      fVmin[2] = 45.5290412902832;
      fVmax[2] = 740.257873535156;
      fVmin[3] = 37.190013885498;
      fVmax[3] = 661.453979492188;
      fVmin[4] = 0.504642367362976;
      fVmax[4] = 3.1229031085968;
      fVmin[5] = 0.870305836200714;
      fVmax[5] = 3.40389847755432;
      fVmin[6] = 50.0762214660645;
      fVmax[6] = 740.257873535156;
      fVmin[7] = 49.8350563049316;
      fVmax[7] = 835.153930664062;

      // initialize input variable types
      fType[0] = 'F';
      fType[1] = 'F';
      fType[2] = 'F';
      fType[3] = 'F';
      fType[4] = 'F';
      fType[5] = 'F';
      fType[6] = 'F';
      fType[7] = 'F';

      // initialize constants
      Initialize();

   }

   // destructor
   virtual ~ReadCFMlpANN() {
      Clear(); // method-specific
   }

   // the classifier response
   // "inputValues" is a vector of input values in the same order as the 
   // variables given to the constructor
   double GetMvaValue( const std::vector<double>& inputValues ) const;

 private:

   // method-specific destructor
   void Clear();

   // common member variables
   const char* fClassName;

   const size_t fNvars;
   size_t GetNvar()           const { return fNvars; }
   char   GetType( int ivar ) const { return fType[ivar]; }

   // normalisation of input variables
   const Bool_t fIsNormalised;
   Bool_t IsNormalised() const { return fIsNormalised; }
   double fVmin[8];
   double fVmax[8];
   double NormVariable( double x, double xmin, double xmax ) const {
      // normalise to output range: [-1, 1]
      return 2*(x - xmin)/(xmax - xmin) - 1.0;
   }

   // type of input variable: 'F' or 'I'
   char   fType[8];

   // initialize internal variables
   void Initialize();
   double GetMvaValue__( const std::vector<double>& inputValues ) const;

   // private members (method specific)
   // not implemented for class: "ReadCFMlpANN"
};
   inline double ReadCFMlpANN::GetMvaValue( const std::vector<double>& inputValues ) const
   {
      // classifier response value
      double retval = 0;

      // classifier response, sanity check first
      if (!IsStatusClean()) {
         std::cout << "Problem in class \"" << fClassName << "\": cannot return classifier response"
                   << " because status is dirty" << std::endl;
         retval = 0;
      }
      else {
         if (IsNormalised()) {
            // normalise variables
            std::vector<double> iV;
            int ivar = 0;
            for (std::vector<double>::const_iterator varIt = inputValues.begin();
                 varIt != inputValues.end(); varIt++, ivar++) {
               iV.push_back(NormVariable( *varIt, fVmin[ivar], fVmax[ivar] ));
            }
            retval = GetMvaValue__( iV );
         }
         else {
            retval = GetMvaValue__( inputValues );
         }
      }

      return retval;
   }
