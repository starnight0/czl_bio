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
// #include <Bpp/Exceptions.h>
// #include <Bpp/Io.all>
// #include <Bpp/Seq.all>
// #include <Bpp/Text.all>
#include <unistd.h>
#ifdef PTHREAD
#include <pthread.h>
#endif

#define CZL_DBG_INFO(c) ( c <<"(" << __FILE__ << ":" << __LINE__ << ":" << __FUNCTION__; )

using namespace std;

namespace czl_bio {

    typedef unsigned char uint8_t;
    typedef unsigned short uint16_t;
    typedef unsigned int uint32_t;
    typedef unsigned long int uint64_t;
    typedef char int8_t;
    typedef short int16_t;
    typedef int int32_t;
    typedef long int int64_t;
    typedef pair<char,string> CSPair;
    typedef pair<string,char> SCPair;

//  typedef int32_t Int;
//  typedef std::pair<Int,Int> Int2;
    
    #define OPT_STR(name, value) (name = value;)
    #define OPT_I(name, value) (name = atol(value.c_str());)
    #define OPT_F(name, value) (name = atof(value.c_str());)
    

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

    void error(const char msg[], const char file[], int line, short is_flush=0);
    void error(string & msg, string & file, int line, short is_flush=0);
    void error(stringstream & msg, const char file[], int line, short is_flush=0);
    void error(stringstream & msg, string & file, int line, short is_flush=0);
    void error(const char msg[], short is_flush=0);
    void error(string & msg, short is_flush=0);
    void error(stringstream & msg_ss, short is_flush=0);

    void warn(const char msg[], const char file[], int line, short is_flush=0);
    void warn(string & msg, string & file, int line, short is_flush=0);
    void warn(stringstream & msg, const char file[], int line, short is_flush=0);
    void warn(stringstream & msg, string & file, int line, short is_flush=0);
    void warn(const char msg[], short is_flush=0);
    void warn(string & msg, short is_flush=0);
    void warn(stringstream & msg_ss, short is_flush=0);

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
    Msg(Msg &);
};
//

/// Initialize a map using ... = CreateMap()()()()...
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

    /* 
     * Option parser
     * @brief represent an option
     * 
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
        string value_v_m; // all values for this option
        bool is_required_m;

        Opt() {}
        Opt(string const & names_str, string const & type, string const & desc, string const & value, const bool is_required)
        {
            StringUtility::str_to_vector(names_str, ',', ' ', name_v);
            type_m = type;
            desc_m = desc;
            value_m = value;
            is_required_m = is_required;
        }

        Opt(vector<string> const & name_v, string const & type, string const & desc, string const & value, const bool is_required)
        {
            name_v_m = name_v;
            type_m   = type;
            desc_m   = desc;
            value_m  = value;
            is_required_m = is_required;
        }

//      Opt(const Opt & src)
//      {
//          name_v_m = src.name_v_m;
//          type_m = src.type_m;
//          desc_m = src.desc_m;
//          is_required_m = src.is_required_m;
//      }

        void print(int indent=2)
        {
            int m = width_m - indent - pre_m -1 ;
            if (m<=0) { indent = 0; }
            char h[indent+pre_m+1];
            int i;
            for (i=0; i<indent; i++) {
                h[i] = ' ';
            }
            h[i] = '\0';
            size_t n = name_v_m.size();
            for (size_t i=1; i<n; i++) {
                cout << h << "--" << std::setw(max_name_ol_m) << std::left << name_v_m[i] << ":\n";
            }
            cout << h << "--" << std::setw(max_name_ol_m) << std::left << name_v_m[0] << " <" << type_m << "> : ";
            if (is_required_m) {
                cout << "*";
            } else {
                cout << " ";
            }
            cout << "\n";
            for (i=indent; i<indent+pre_m; i++) h[i] = ' ';
            h[i] = '\0';
            int m0 = 0, m1 = m;
            while (m0<desc_m.size()) {
                cout << h << desc_m.substr(m0, m) << "\n";
                m0 += m;
            }
        }

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
    
    /*
     * Parsing Options
     * @brief
     *
     * every option has a name, type, description, required flag, list of values
     * variabe name is the same as the option name
     * 
     * auto read option from configure files and command line
     * 
     * @date 2015-10-27
     *
     */
    int get_opts( vector<Opts> opt_v, int argc, char** argv )
    {
        map<string, int> opt_sl;
        for (int i=0; i<opt_v.size(); i++) {
            Opt & p = opt_v[i];
            for (int j=0; j<p.name_v_m.size(); j++) {
                if ( opt_sl.find(p.name_v_m[j])!=opt_sl.end() ) {
                    cerr << "E: Duplicate option: '" << p.name_v_m[j] << "'";
                    CZL_DBG_INFO(cerr);
                    return 1;
                }
                opt_sl[ p.name_v_m[i] ] = i;
            }
        }
        /// get option
        // {{{
        for (k=1; k<argc; k++) {
            string name;
            string value;
            if ( argv[k][0]=='-' ) {
                if ( argv[k][1]=='-' ) {
                    name = argv[k]+2;
                } else {
                    name = argv[k]+1;
                }
                if ( name=="h" || name == "help") {
                    print_usage();
                    #ifdef PTHREAD
                    pthread_exit(NULL);
                    #else
                    return 0;
                    #endif
                } else if ( name=="v" || name == "version") {
                    cout << prog_version << endl;
                    #ifdef PTHREAD
                    pthread_exit(NULL);
                    #else
                    return 0;
                    #endif
                }

                for (i=0; i<name.size(); i++) {
                    if (name[i]=='-') name[i]='_';
                }
                if ( opt_sl.find(name)!=opt_sl.end() ) name = opt_sl[name];
                if ( opt_n.find(name)==opt_n.end() ) {
                    cerr << "Error: Unavailable option '" << name << "'" << endl;
                    #ifdef PTHREAD
                    pthread_exit(NULL);
                    #else
                    return 0;
                    #endif
                } else {
                    if (opt_n[name]==1) {
                        k++;
                        if (k>=argc) {
                            cerr << "Error: option " << name << "need value." << endl;
                            #ifdef PTHREAD
                            pthread_exit(NULL);
                            #else
                            return 0;
                            #endif
                        }
                        value.assign(argv[k]);
                    }
                }
                opt[name] = value;
            } else {
                cerr << "Error: option not start with '-': " << argv[k] << endl;
                #ifdef PTHREAD
                pthread_exit(NULL);
                #else
                return 0;
                #endif
            }
        }
        /// }}}

        if ( opt.find("config")!=opt.end() ) { // -c OR --config
            /// read configure file
            /// {{{
            conf_file = opt["config"];
            opt.erase("config");
            cerr << "Read from configure file " << conf_file << endl;
            if (!conf_file.empty()) {
                fs.open(conf_file.c_str(), fstream::in);
                if (fs.fail()) {
                    cerr << "Can't read configure file " << conf_file << endl;
                    return 1;
                }
                while (!fs.eof()) {
                    getline(fs, str);
                    boost::trim(str);

                    if (str.empty()) { continue; }

                    if (str[0]=='[') { continue; }
                    if (str[0]=='#') { continue; }

                    i=str.find_first_of('=');
                    if (i == std::string::npos) continue;
                    
                    string name, value;
                    name=str.substr(0, i);
                    value=str.substr(i+1);

                    if (name.empty()) continue;
                    boost::trim(name);
                    if (value.empty()) continue;
                    boost::trim(value);
                    if ( opt_sl.find(name)!=opt_sl.end() ) name = opt_sl[name];
                    if ( opt_n.find(name)!=opt_n.end() ) {
                        opt[name] = value;
                    } else {
                        cerr << "Can't recognize option " << name << " = " << value << " Skip." << endl;
                    }
                }
                fs.close();
            }
            /// }}}
        }
    }

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


    /* String Utility
     */
    class StringUtility {
    public:
        static int get_next(const string & line, const char sep, const int pos, string sep, string & out);
        static int skip(string & line, int pos, string sep, int n=1);
        static int trim_left(string & str, const char trim);
        static int trim_right(string & str, const char trim);
        static int trim(string & str, const char trim);
        static int str_to_vector(const string & str, const char sep, vector<string> & out_v);
        static int str_to_vector(const string & str, const char sep, const char trim, vector<string> & out_v);
    };

};

#endif
