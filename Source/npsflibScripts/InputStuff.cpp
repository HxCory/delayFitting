#include <string.h>
#include <stdlib.h>
#include "InputStuff.h"


extern Parameter<string> outputFolder; // it's defined outside

ConfigFileParser cfg; // global


//void create_dir(std::string & name){
inline std::string create_dir(std::string name){
  /**
   * name - name of the directory to be created
   *        (for instance a global parameter)
   * In case directory name exists the function 
   * tries name__1, name__2, ... etc. 
   * 
   * If it was succesfull returns the new directory name
   * ( and also used to append it to the file "_jobnames" )
   *
   * Exits the program in case of failure.
   * 
   * possible simple modifications:
   * - change cerr to other stream, like log file
   *   (now both success and different failures are reported),
       change format of logs
   * - change the rule for choosing names
   * - change exit to other action
   * - change permissions of created directory 
   * - set if you want to respect umask or not (set const bool parameter) 
   * - maximal nr of names it tries is 10000 now
   */
  using namespace std;
  int i=0;
  mode_t mask;
  FILE *jobnames;
  const bool respect_umask=false;
  stringstream tmp;
  
  tmp.str(""); // set string to empty
  tmp<<name;
  if(!respect_umask) mask=umask(007); //optionally skip system-wide umask
  while(i<10000){
    // here permissions, see sys/stat.h for names (may be google)
    if(mkdir( tmp.str().c_str(),S_IRWXU|S_IRGRP|S_IXGRP|S_IROTH|S_IXOTH)==0){
      //after successful creation of directory do the following stuff 
      //and exit from this function returning the new name
      /*
      jobnames=fopen("_jobnames","a");
      if(jobnames==NULL){
	cerr<<__FILE__<<":"<<__LINE__<<": "<<strerror(errno)<<"; _jobnames not updated"<<endl;
      }
      else{
	fprintf(jobnames,"%s\n",tmp.str().c_str());
	fclose(jobnames);
      }
      */
      //here cerr
      cout<<__FILE__<<":"<<__LINE__<<": created directory \""<<tmp.str()<<"\""<<endl;
      //if(!respect_umask) umask(mask); // reset umask  at the return
      return tmp.str();
    }
    if(errno!=EEXIST) break;
    i++;
    tmp.str(""); //reset the string
    // here creation of suffixes
    tmp<<name<<"__"<<i;
  }
  //the following is in case of failure
  cerr<<__FILE__<<":"<<__LINE__<<": "<<strerror(errno)<<". Were not able to create directory "<< name<<endl;
  cerr<<"exit"<<endl;
  //if(!respect_umask) umask(mask); // reset umask  at the return
  exit(1);
}

void touchdir(std::string name){
  /**
   * name - name of the directory to be created
   *        (for instance a global parameter)
   *
   * Exits the program in case of failure.
   * 
   */
  using namespace std;
  mode_t mask;
  const bool respect_umask=false;

  // check exsistence of dir
  struct stat st;
  stat(name.c_str(), &st);
  if(S_ISDIR(st.st_mode)) return;  // directory already there
  
  if(!respect_umask) mask=umask(0); //optionally skip system-wide umask
  // here permissions, see sys/stat.h for names (may be google)
  if(mkdir( name.c_str(),S_IRWXU|S_IWGRP|S_IRGRP|S_IXGRP|S_IROTH|S_IXOTH)==0){
    //here cerr
    cout<<__FILE__<<":"<<__LINE__<<": created directory \""<<name<<"\""<<endl;
    if(!respect_umask) umask(mask); // reset umask  at the return
  } else {
    cerr<<__FILE__<<":"<<__LINE__<<": "<<strerror(errno)<<". Were not able to create directory "<< name<<endl;
    if(!respect_umask) umask(mask); // reset umask  at the return
    exit(1);
  }
}

void parseInput(int argc, char *argv[]){
  string cfgName("example.cfg");
  if( argc == 2 ) {
    cfgName=argv[1];
  } else {
    cerr<<"usage: "<<argv[0]<<" configfile"<<endl;
    cerr<<"trying default: "<<cfgName<<endl;
  }
  /***************************************/
  // read configuration file, name may be obtained from command line
  cfg.read( cfgName.c_str() );
  // the next line gets a pointer to the ParameterMap object:
  //ParameterMap& parameters = ParameterMap::instance();
  ParameterMap parameters;
  // initialize all Paramter<T> type objects from the config file
  parameters.init( cfg );
  // | could also simply become global::paramters.init( cfg )    |
  // | instead of the two lines above                            |

  // create the output directory (after reading the name of the job from cfg)
  string job = parameters["jobname"]->getValueString();
  outputFolder.set(create_dir( job )+"/");      // with slash at the end!

  cout<<parameters;
  cout<<"Press Ctrl-C to stop"<<endl;
  //sleep(3);

  /****************************************/
}



