#ifndef CZL_COMMON_H
#define CZL_COMMON_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cctype>
#include <cfloat>
#include <cmath>
#include <cstdlib>
#include <cassert>
#include <vector>
#include <queue>
#include <list>
#include <string>
#include <map>
#include <algorithm>
#include <bitset>
#include <set>
#include <utility>
#include <iterator>
#include <iomanip>
#include <exception>
#include <stdexcept>
#include <sys/stat.h>
#include <cerrno>
// #include <Bpp/Exceptions.h>
// #include <Bpp/Io.all>
// #include <Bpp/Seq.all>
// #include <Bpp/Text.all>
#include <unistd.h>
#ifdef PTHREAD
#include <pthread.h>
#endif

//#define CZL_DBG_INFO "(" << __FILE__ << ":" << __LINE__ << ":" << __FUNCTION__ << ")"
#define CZL_DBG_INFO (string("(") + __FILE__ + ":" + itos(__LINE__) + ":" + __FUNCTION__ + ")")
#define CZL_CUR_TIME ( str_time(time(NULL)) + ", " + ftos((float)clock()/CLOCKS_PER_SEC ) )

using namespace std;

namespace czl_bio {

#ifndef _STDINT_H
    typedef unsigned char      uint8_t;
    typedef unsigned short     uint16_t;
    typedef unsigned int       uint32_t;
    typedef unsigned long int  uint64_t;
    typedef char  int8_t;
    typedef short int16_t;
    typedef int   int32_t;
#if __WORDSIZE == 64
    typedef long int      int64_t;
#else
    typedef long long int int64_t;
#endif
#endif
    typedef pair<char,string> CSPair;
    typedef pair<string,char> SCPair;

    class StringUtility;

//  typedef int32_t Int;
//  typedef std::pair<Int,Int> Int2;

    #define DECL_TYPE_NAME(x) template<> struct type_name<x> { static const char* name() {return #x;} }   

	string atos(char *str);
    char* itoa(int n);
    string itos(int n);
    char* ltoa(long int n);
    string ltos(long int n);
    string ftos(double f);
    string ftos(float f);
    char* czl_itoa(int i, int base=10);
    string czl_itos(int i, int base=10);
    char* czl_ltoa(long int i, int base=10);
    string czl_ltos(long int i, int base=10);
    char *czl_ftoa(double f, int digit, int base=10);
    string czl_ftos(double f, int digit, int base=10);

    /// range union, intersect
    // {{{
    template<typename T>
    int czl_union(T & begin1, T & end1, T & begin2, T & end2, T & out_begin, T & out_end)
    {
        if ( begin1 < begin2 ) {
            out_begin = begin1;
        } else {
            out_begin = begin2;
        }
        if ( end1 < end2 ) {
            out_end = end2;
        } else {
            out_end = end1;
        }
        return 0;
    }

    template<typename T>
    int czl_union_len(T & begin1, T & end1, T & begin2, T & end2)
    {
        T b, e;
        czl_union(begin1, end1, begin2, end2, b, e);
        return e-b;
    }

    template<typename T>
    int czl_intersect(T & begin1, T & end1, T & begin2, T & end2, T & out_begin, T & out_end)
    {
        T b, e;
        if ( begin1 < begin2 ) {
            out_begin = begin2;
        } else {
            out_begin = begin1;
        }
        if ( end1 < end2 ) {
            out_end = end1;
        } else {
            out_end = end2;
        }
        if ( out_end < out_begin ) { out_begin=0; out_end=0; }
        return 0;
    }

    template<typename T>
    int czl_intersect_len(T & begin1, T & end1, T & begin2, T & end2)
    {
        T b, e;
        czl_intersect(begin1, end1, begin2, end2, b, e);
        if ( e-b>0 ) return e-b;
        else return 0;
    }
    // }}}

    /* 
     * @brief Clasee Msg: Output Error, Warn Message
     */
    class Msg {
    protected:
    //    static const uint32_t N;
    //    static uint32_t nw; // number of warn
    //    static int32_t m;
    //    static std::vector<std::string> msg_m;
        string file_m;
        ofstream os_m;
    //    short is_flush;
        short is_init_m;
#ifdef PTHREAD
        pthread_mutex_t lock_m;
#endif
    public:
        Msg();
        Msg(string & file);
        Msg(const char file[]);
        ~Msg();
        int init(string & file);
        int init(const char file[]);
        short is_init();
        short is_good();
        short good();
        short is_fail();
        short fail();

		Msg & print(const char head[], const char msg[], short is_flush=0, short where=3);
		Msg & print(const string & head, const string & msg, short is_flush=0, short where=3);
		Msg & print(const char msg[], short is_flush=0, short where=3);
		Msg & print(const string & msg, short is_flush=0, short where=3);

        Msg & error(const char msg[], short is_flush=0, short where=3);
        Msg & error(const char msg[], const char file[], int line, short is_flush=0, short where=3);
        Msg & error(string const & msg, string & file, int line, short is_flush=0, short where=3);
        Msg & error(stringstream & msg, const char file[], int line, short is_flush=0, short where=3);
        Msg & error(stringstream & msg, string & file, int line, short is_flush=0, short where=3);
        Msg & error(string const & msg, short is_flush=0, short where=3);
        Msg & error(stringstream & msg_ss, short is_flush=0, short where=3);

        Msg & warn(const char msg[], short is_flush=0, short where=3);
        Msg & warn(const char msg[], const char file[], int line, short is_flush=0, short where=3);
        Msg & warn(string const & msg, string & file, int line, short is_flush=0, short where=3);
        Msg & warn(stringstream & msg, const char file[], int line, short is_flush=0, short where=3);
        Msg & warn(stringstream & msg, string & file, int line, short is_flush=0, short where=3);
        Msg & warn(string const & msg, short is_flush=0, short where=3);
        Msg & warn(stringstream & msg_ss, short is_flush=0, short where=3);

        Msg& operator << (const int i);
        Msg& operator << (const unsigned int i);
        Msg& operator << (const long int i);
        Msg& operator << (const unsigned long int i);
        Msg& operator << (const char c);
        Msg& operator << (const float f);
        Msg& operator << (const double f);
        Msg& operator << (const char *str);
        Msg& operator << (const string & s);
        Msg& operator << (streambuf* sb );
        Msg& operator << (ostream& (*pf)(ostream&));
        Msg& operator << (ios& (*pf)(ios&));
        Msg& operator << (ios_base& (*pf)(ios_base&));

        Msg& endl(Msg & msg);
        Msg& flush();
        Msg& flush(Msg& msg);
    private:
        Msg(Msg const &);
    };
	static Msg msg;
    //

    /*
     * class CreateMap
     *
     * For initializing a map outside main()
     * Initialize a map using ... = CreateMap()()()()...
     */
    // {{{
    template <typename T, typename U>
    class CreateMap
    {
    private:
        std::map<T, U> map_m;
    public:
        CreateMap(const T& key, const U& val)
        {
            map_m[key] = val;
        }

        CreateMap<T, U>& operator()(const T& key, const U& val)
        {
            map_m[key] = val;
            return *this;
        }

        operator std::map<T, U>()
        {
            return map_m;
        }
    };
    // }}}

    /*
     * class Range
     */
    template <typename Pos_T>
    class Range: public pair<Pos_T,Pos_T> {
    public:
        Range(): pair<Pos_T,Pos_T>(), strand_m('?') {}
        Range(const Range & src): pair<Pos_T, Pos_T>(src)
        {
            strand_m = src.strand_m;
        }

        Range & operator = (const Range & src)
        {
            pair<Pos_T, Pos_T>::operator=(pair<Pos_T,Pos_T>(src));
            strand_m = src.strand_m;
            return *this;
        }

        Pos_T get_len() {
            return this->second - this->first;
        }
        Pos_T get_begin() {
            return this->first;
        }
        Pos_T get_end() {
            return this->second;
        }
        char get_strand() {
            return strand_m;
        }

        void set_begin(Pos_T begin) {
            this->first = begin;
        }
        void set_end(Pos_T end) {
            this->second = end;
        }
        void set_strand(char strand) {
            strand_m = strand;
        }
        void set(Pos_T begin, Pos_T end, char strand)
        {
            this->first    = begin;
            this->second   = end;
            strand_m = strand;
        }

        void clear()
        {
            memset(this, 0, sizeof(*this));
            strand_m = '?';
        }

        ~Range()
        {
        }
    private:
        char strand_m;
    };

    /**
     * @brief String Utility
     */
    class StringUtility {
    public:
        static int find(const string & line, const char sep, const int pos, string & out);
        static int find(const string & line, string const & sep, const int pos, string & out);
        static int find_first_of(const string & line, string const & sep, const int pos, string & out);
        static int skip(string & line, int pos, string sep, int n=1);
        static int trim_left(string & str, const char trim);
        static int trim_left(string & str, string const & trim);
        static int trim_right(string & str, const char trim);
        static int trim_right(string & str, string const & trim);
        static int trim(string & str, const char trim);
        static int trim(string & str, string const & trim);
        static int split(const string & str, const char sep, vector<string> & out_v);
        static int split(const string & str, const string & sep, vector<string> & out_v);
        static int split(const string & str, const char sep, const char trim, vector<string> & out_v);
        static int split(const string & str, string const & sep, const char trim, vector<string> & out_v);
        static int split_first_of(const string & str, string const & sep, vector<string> & out_v);
        static int split_first_of(const string & str, string const & sep, const char trim, vector<string> & out_v);
        static int erase_all(string &str, const string & set);

        static void print(const string & str, int indent=0, int width=0, ostream & os=std::cout);
        static void to_upper(string & str);
        static void to_upper_copy(string const & str, string & out);
        static void to_lower(string & str);
        static void to_lower_copy(string const & str, string & out);

        // format printing
        // {{{
        /**
         * @brief  format printing from one string to another string
         * @param  s   input string
         * @param  w   window width
         * @param  lm  left margin
         * @param  rm  right margin
         * @param  out output string
         * @return  ref to the output string
         * @remark  The output string is not clean first,
         *          output string is appended. 
         * @sa  ostream & print(string const &, int, int, int, ostream &)
         */
        static string & mar_print(string const &s, int w, int lm, int rm, string & out);

        /**
         * @brief  format printing from one string to output stream
         * @param  s   input string
         * @param  w   window width
         * @param  lm  left margin
         * @param  rm  right margin
         * @param  out output stream
         * @return  ref to the output stream
         * @sa  string & print(string const &, int, int, int, string &)
         */
        static ostream & mar_print(string const &s, int w, int lm, int rm, ostream & out);
        // }}}
    };

    /**
     * @brief represent an option
     */
    // {{{
    class Opt {
    public:
        static const int max_name_ol_m = 20;
        static const int width_m = 80;
        static const int pre_m = 2;

        vector<string> name_v_m;
        string type_m;
        string desc_m;
        string value_m; // all values for this option
        int n_m;  // number of parameter
        bool is_required_m;
        bool is_set_m;

        /**
         * @brief Constructor using C-type chars
         * @param  names_str  name of the option (don't include hyphen at the 
         *                    head) 
         * @param  desc       description of the option (will show in usage 
         *                    output)
         * @param  value      unparsed values of the option, use for set 
         *                    default value
         * @param  n          number of parameter for the option
         * @param  is_required  is the option required for the program
         */
        Opt(char const * names_str, char const * type, char const * desc, char const * value, int n, bool is_required);

        /**
         * @brief Constructor using C++ string
         * @param  names_str  name of the option (don't include hyphen at the 
         *                    head) 
         * @param  desc       description of the option (will show in usage 
         *                    output)
         * @param  value      unparsed values of the option, use for set 
         *                    default value
         * @param  n          number of parameter for the option
         * @param  is_required  is the option required for the program
         * @sa  Opt(char const *, char const *, char const *, char const *, 
         *      int, bool)
         */
        Opt(string const & names_str, string const & type, string const & desc,
                string const & value, int n, bool is_required);

        /**
         * @brief Constructor using C++ string
         * @param  name_v  vector of names of the option (don't include hyphen
         *                 at the head) 
         * @param  desc    description of the option (will show in usage 
         *                 output)
         * @param  value   unparsed values of the option, use for set 
         *                 default value
         * @param  n       number of parameter for the option
         * @param  is_required  is the option required for the program
         * @sa     Opt(char const *, char const *, char const *, char const *,
         *         int, bool), Opt(string const, string const &, 
         *         string const & , string const &, int, bool)
         */
        Opt(vector<string> const & name_v, string const & type, 
            string const & desc, string const & value, int n, 
            const bool is_required);

//      Opt(const Opt & src)
//      {
//          name_v_m = src.name_v_m;
//          type_m = src.type_m;
//          desc_m = src.desc_m;
//          is_required_m = src.is_required_m;
//      }

        virtual int str_to_value() { return 0; }
    //  virtual int value_to_str() {};

        /**
         * @brief print the option usage
         * @param  os  Output stream
         * @param  indent  Indent length
         */
        void print_usage(ostream & os, int indent=2);

        /**
         * @brief print the option in format: NAME = VALUE
         * @param  os  Output stream
         * @param  indent  Indent length
         */
        void print_config(ostream & os, int indent=0);

        /**
         * @brief  check is there a value for the option, return true if the 
         *         the value is not empty
         * @return  true if the value is not empty
         */
        bool has_value()
        {
            return !value_m.empty();
        }

        /**
         * @brief  check if the option is set from command line or config file
         * @return  true if the option is set from command line or config file.
         */
        bool is_set()
        {
            return is_set_m;
        }

        virtual ~Opt () {}
    private:
        Opt & operator=(const Opt & src) /// copy is not available
        {
            name_v_m = src.name_v_m;
            type_m = src.type_m;
            desc_m = src.desc_m;
            is_required_m = src.is_required_m;
            return *this;
        }

    };
    // }}}

    /**
     * @brief Represent command line option
     * @tparam  T  the type of the option that will be changed to
     * @tparam  str_to_value_func  function to convert a string to type T
     */
    template<typename T, int (*str_to_value_func)(const string &, T &)>
    class OptT: public Opt {
    private:
        T v_m;
    public:
        /**
         * @brief Constructor using c-type char
         * @param  names_str  name of the option (don't include hyphen at the 
         *                    head) 
         * @param  desc       description of the option (will show in usage 
         *                    output)
         * @param  value      unparsed values of the option, use for set 
         *                    default value
         * @param  n          number of parameter for the option
         * @param  is_required  is the option required for the program
         * @sa  OptT(string const &, string const &, string const &, 
                string const &, int, bool)
         * @date 2015-11-30
         */
        OptT(char const * names_str, char const *type, char const * desc,
                char const * value, int n, bool is_required): 
                Opt(names_str, type, desc, value, n, is_required)
        {
            if (!value_m.empty()) {
                str_to_value();
            }
        }

        /**
         * @brief Constructor using string
         * @param  names_str  name of the option (don't include hyphen at the 
         *                    head) 
         * @param  desc       description of the option (will show in usage 
         *                    output)
         * @param  value      unparsed values of the option, use for set 
         *                    default value
         * @param  n          number of parameter for the option
         * @param  is_required  is the option required for the program
         * @sa  OptT(char const *, char const *, char const *, char const *,
                int, bool)
         * @date 2015-11-30
         */
        OptT(string const & names_str, string const & type, string const & desc, 
                string const & value, int n, bool is_required):
                Opt(names_str, type, desc, value, n, is_required)
        {
            if (!value_m.empty()) {
                str_to_value();
            }
        }


        int str_to_value()
        {
            if (str_to_value_func==NULL) return 0;
            return (*str_to_value_func)(value_m, v_m);
        }

        /**
         * @brief  get the converted value
         * @return  The refence to the converted value;
         */
        T & get_value_ref()
        {
            return v_m;
        }

        ~OptT() {}
    };

    /**
     * @brief Parsing a group of OptT-type options
     *
     * every option has a name, type, description, required flag, list of values
     * variabe name is the same as the option name
     * 
     * auto read option from configure files and command line
     * 
     * @date 2015-10-27
     *
     */
    class Opts {
    private:
        int argc_m;
        char** argv_m;
        vector<Opt*> opt_v_m;
        map<string, int> opt_sl_m;

    public:
        /**
         * @brief  Constructor. It's create two default options
         *         'help' , 'version' and 'config'
         */
        Opts();

        /**
         * @brief  Create an option and append, may throw and catch bad_alloc 
         *         exception if can't allocate the option
         * @tparam  OT  Any type of OptT, derived from Opt
         * @param  names_str  name of the option (don't include hyphen at the 
         *                    head) 
         * @param  desc       description of the option (will show in usage 
         *                    output)
         * @param  value      unparsed values of the option, use for set 
         *                    default value
         * @param  n          number of parameter for the option
         * @param  is_required  is the option required for the program
         * @return  0 if success, 1 if there is duplicated option names,
         *          may throw and catch bad_alloc exception
         * @remark  Don't DELETE or FREE the returned option
         * @see     OptT::OptT(string const, string const, string const, int, bool)
         */
        template<typename OT>
        OT* create(string const & names_str, string const & type, string const & desc, string const & value, int n, bool is_required)
        {
            OT *p = new OT(names_str, type, desc, value, n, is_required);
            if (p!=NULL) {
                if ( push_back(p) ) { // fail to push
                    delete p;
                    p = NULL;
                }
            }
            return p;
        }

        /**
         * @brief   get an option by name
         * @tparam  OT    any OptT class
         * @param   name  name of the option
         * @return  Pointer to the option if the option name exists;
         *          NULL otherwise.
         */
        template<typename OT>
        OT* get_opt_by_name(const string & name)
        {
            if ( opt_sl_m.find(name)==opt_sl_m.end() ) {
                return NULL;
            } else {
                return dynamic_cast<OT*>(opt_v_m[ opt_sl_m[name] ]);
            }
        }


        /**
         * @brief  Create an option and append, may throw and catch bad_alloc 
         *         exception if can't allocate the option
         * @tparam  OT  Any type of OptT, derived from Opt
         * @param  names_str  name of the option (don't include hyphen at the 
         *                    head) 
         * @param  desc       description of the option (will show in usage 
         *                    output)
         * @param  value      unparsed values of the option, use for set 
         *                    default value
         * @param  n          number of parameter for the option
         * @param  is_required  is the option required for the program
         * @return  0 if success, 1 if there is duplicated option names,
         *          may throw and catch bad_alloc exception
         * @remark  Don't DELETE or FREE the returned option
         * @see     OptT::OptT(char const *, char const * , char const * , int , bool )
         */
        template<typename OT>
        OT* create(char const * names_str, char const * desc, char const * value, int n, bool is_required)
        {
            OT *p = new OT(names_str, desc, value, n, is_required);
            if (p!=NULL) {
                if ( push_back(p) ) { // fail to push
                    delete p;
                    p = NULL;
                }
            }
            return p;
        }


        /**
         * @brief parse options from command line
         *
         * Parse a shell command line 
         * 'PROGRAM --name1 value1 --name2 value2 and so on'.
         * 
         * Command line can contain a config file use '--config' or '-c'
         * option, which will auto parse. Option name can occur more than 
         * once, but only the last one will be used except the '-c' option,
         * All config file will be parsed, and the later value will cover the 
         * former value.
         * 
         * You'll get a warning when use same name more than once,
         *
         * @param  argc  number of options in command line, include the 
         *               program name, same as 'argc' in C/C++ main
         * @param  argv  raw options in command line, include the 
         *               program name, same as 'argv' in C/C++ main
         * @return  0 if success, 
         *          1 if the option in command line is not avaliable
         *          2 if the config file is not available
         *          3 if the option in config file is not avaliable
         *          4 if the option don't have value in command line;
         *            but empty name and value in config file will be ignore
         */
        int parse(int argc, char** argv);

        /**
         * @brief check if all required option set
         * @return  true if all required otpions are set, false otherwise.
         */
        bool is_all_set();
        /**
         * @brief  print the usage of all options
         * @param  s  Output stream to print
         */
        void print_usage(ostream & s);

        /**
         * @brief  print the option NAME = VALUE pair
         * @param  s  Output stream to print
         */
        void print_config(ostream & s);

        static int str_to_value_b(string const & str, bool & v) 
        {
            if (str.size()<=0) return -1;
            else {
                string str1 = str;
                std::transform(str.begin(), str.end(), str1.begin(), ::tolower);
                if (str1=="t" || str1=="true") v=true;
                else v=false;
            }
            return 0;
        }

        static int str_to_value_c(string const & str, char & v) 
        { if (str.size()<=0) return -1; else v = str[0]; return 0;}

        static int str_to_value_i(string const & str, int & v) 
        { v = atoi(str.c_str()); return 0;}

        static int str_to_value_l(string const & str, long int & v) 
        { v = atol(str.c_str()); return 0;}

        static int str_to_value_f(string const & str, float & v) 
        { v = atof(str.c_str()); return 0;}

        static int str_to_value_lf(string const & str, double & v) 
        { v = atof(str.c_str()); return 0;}

        static int str_to_value_s(string const & str, string & v) 
        { v.assign(str); return 0;}

        static int str_to_value_iv(string const & str, vector<int> & v);
        static int str_to_value_lv(string const & str, vector<long> & v);
        static int str_to_value_fv(string const & str, vector<float> & v);
        static int str_to_value_lfv(string const & str, vector<double> & v);
        static int str_to_value_sv(string const & str, vector<string> & v);

        /**
         * Destructor
         */
        ~Opts();
    private:
        /**
         * @brief  append an option, may throw and catch bad_alloc exception
         *         users are not allowed to push an option using this function,
         *         Use 'create' instea.
         * @param  p  pointer to the option
         * @return  0 if success, 1 if there is duplicated option names,
         *          may throw and catch bad_alloc exception
         * @see    create()
         */
        int push_back(Opt* p);
    };

    /// bool Option
    typedef OptT<bool, Opts::str_to_value_b> OptB;
    /// char Option
    typedef OptT<char, Opts::str_to_value_c> OptC;
    /// Integer Option
    typedef OptT<int, Opts::str_to_value_i> OptI;
    typedef OptT<long, Opts::str_to_value_l> OptL;
    /// Float/Double Option
    typedef OptT<float, Opts::str_to_value_f> OptF;
    typedef OptT<double, Opts::str_to_value_lf> OptLF;
    /// String Option
    typedef OptT<std::string, Opts::str_to_value_s> OptS;
    /// Integer Array Option
    typedef OptT<vector<int>, Opts::str_to_value_iv> OptIV;
    typedef OptT<vector<long>, Opts::str_to_value_lv> OptLV;
    /// Double Array Option
    typedef OptT<vector<float>, Opts::str_to_value_fv> OptFV;
    typedef OptT<vector<double>, Opts::str_to_value_lfv> OptLFV;
    /// String Array Option
    typedef OptT<vector<std::string>, Opts::str_to_value_sv> OptSV;


    /*
    void print_opts( vector<Opt> & opt_v )
    {
        for (vector<Opt>::iterator it=opt_v.begin(); it!=opt_v.end(); it++) {
            it->print();
        }
    }

    void print_opts_value( ostream & s, vector<Opt> & opt_v )
    {
        for (vector<Opt>::iterator it=opt_v.begin(); it!=opt_v.end(); it++) {
            s << std::setw(Opt::max_name_ol_m) << std::left << it->name_v_m[0] << "= ";
            for (int i=0; i<it->value_v_m.size(); i++) {
                s << it->value_v_m[i] << ";";
            }
            s << "\n";
        }
        s << std::flush;
    }
    */

    /**
     * @brief Some file operating function
     */
    class File {
    public:
        /**
         * @brief  Check if a regular file exists
         */
        static bool exists_file(string const &fname);
        /**
         * @brief  Check if a directory exists
         */
        static bool exists_dir(string const &fname);
        /**
         * @brief  Create a new directory
         */
        static int create_dir(string const &fname);
    };

    /**
     * @brief  check if big-endian
     * @return  true if the machine is big-endian, false otherwise
     */
    bool is_be();

    /**
     * @brief  get the number of address bits
     * @return  the number of address bits on the current machine
     */
    int machine_bit();

    /**
     * @brief print local time to a string
     */
    string str_time(time_t const & t);

    /**
     * @brief print difftime to a string (day,hour,minute,second)
     */
    string str_difftime(time_t const & end, time_t const & begin);

template<typename T>
void print_array(int n1, int n2, T** a)
{
    for (int i1=0; i1<n1; i1++) {
        for (int i2=0; i2<n2; i2++) {
            cout << a[i1][i2] << "  ";
        }
        cout << "\n";
    }
}

template<typename T>
void print_array(vector< vector<T> > const & a)
{
    int n1 = a.size();
    for (int i1=0; i1<n1; i1++) {
        for (int i2=0; i2<a[i1].size(); i2++) {
            cout << a[i1][i2] << "  ";
        }
        cout << "\n";
    }
}

template<typename T>
void print_array(ostream & s, int n1, int n2, T** a)
{
    for (int i1=0; i1<n1; i1++) {
        for (int i2=0; i2<n2; i2++) {
            s << a[i1][i2] << "  ";
        }
        s << "\n";
    }
}

template<typename T>
void print_array(ostream & s, vector< vector<T> > const & a)
{
    int n1 = a.size();
    for (int i1=0; i1<n1; i1++) {
        for (int i2=0; i2<a[i1].size(); i2++) {
            s << a[i1][i2] << "  ";
        }
        s << "\n";
    }
}

};

#endif
