/****************************************************************************** 
 * @file  czl_common.cpp
 * @brief Common Functions
 * @author ZELIN CHEN
 ******************************************************************************/
#include "czl_common.hpp"
#include <ctype.h>

using namespace std;

namespace czl_bio {

/*
 * string -- number convertor
 */
// {{{
string atos(char *str)
{
    return string(str);
}

char* itoa(int n)
{
    int base=10;
    int sign=n;
    int l=0;
    int m=n;
    if (sign < 0) n = -n;
    do {
        m=m/base;
        l++;
    } while (m>0);
    char *s = NULL;
    if ( sign>=0 ) {
        s = new char[l+1];
    } else {
        l++;
        s = new char[l+1];
    }
    sprintf(s, "%i", n);
    return s;
}

string itos(int n)
{
    stringstream ss;
    ss << n;
    return ss.str();
}

char* ltoa(long int n)
{
    long int sign=n;
    int l=0;
    long int m=n;
    if (sign < 0) n = -n;
    do {
        m=m/10;
        l++;
    } while (m>0);
    char *s = NULL;
    if ( sign>=0 ) {
        s = new char[l+1];
    } else {
        l++;
        s = new char[l+1];
    }
    sprintf(s, "%ld", n);
    return s;
}

string ltos(long int n)
{
    stringstream ss;
    ss << n;
    return ss.str();
}

string ftos(double f)
{
    stringstream ss;
    ss << f;
    return ss.str();
}

string ftos(float f)
{
    stringstream ss;
    ss << f;
    return ss.str();
}


/// @brief change int to C char
char* czl_itoa(int n, int base)
{
    int sign=n;
    int l=0;
    int m=n;

    if (base<36) return NULL;

    if (sign < 0) n = -n;
    do {
        m=m/base;
        l++;
    } while (m>0);
    char *s = NULL;
    if ( sign>=0 ) {
        s = new char[l+1];
    } else {
        l++;
        s = new char[l+1];
    }
    s[l--] = '\0';
    do {       // generate digits in reverse order
        m = n % base;
        if (base >= 10 && m >= 10) {
            s[l--] = m-10 + 'A';
        } else {
            s[l--] = m + '0';
        }
        n /= base;
    } while (n > 0);
    if (sign < 0) s[l--] = '-';

    return s;
}

string czl_itos(int n, int base)
{
    char * s=czl_itoa(n, base);
    string str;
    if (s!=NULL) {
        str.assign(s);
        delete s;
    }
    return str;
}

char* czl_ltoa(long int n, int base)
{
    long int sign=n;
    int l=0;
    long int m=n;

    if (base<36) return NULL;

    if (sign < 0) n = -n;
    do {
        m=m/base;
        l++;
    } while (m>0);
    char *s = NULL;
    if ( sign>=0 ) {
        s = new char[l+1];
    } else {
        l++;
        s = new char[l+1];
    }
    s[l--] = '\0';
    do {       // generate digits in reverse order
        m = n % base;
        if (base >= 10 && m >= 10) {
            s[l--] = m-10 + 'A';
        } else {
            s[l--] = m + '0';
        }
        n /= base;
    } while (n > 0);
    if (sign < 0) s[l--] = '-';

    return s;
}

string czl_ltos(long int n, int base)
{
    char * s=czl_ltoa(n, base);
    string str;
    if (s!=NULL) {
        str.assign(s);
        delete s;
    }
    return str;
}

char *czl_ftoa(double f, int digit, int base)
{
    long int n=0;
    for (int i=0; i<digit; i++) {
        f *= base;
    }
    n = static_cast<long int>(f);

    long int sign=n;
    int l=0;
    long int m=n;

    if (base<36) return NULL;

    if (sign < 0) n = -n;
    do {
        m=m/base;
        l++;
    } while (m>0);
    char *s = NULL;
    if ( digit>0 ) l++;
    if ( sign>=0 ) {
        s = new char[l+1];
    } else {
        l++;
        s = new char[l+1];
    }
    s[l--] = '\0';
    do {       // generate digits in reverse order
        m = n % base;
        if (base >= 10 && m >= 10) {
            s[l--] = m-10 + 'A';
        } else {
            s[l--] = m + '0';
        }
        if (digit>0) {
            digit--;
            if (digit==0) {
                s[l--]='.';
            }
        }
        n /= base;
    } while (n > 0);
    if (sign < 0) s[l--] = '-';

    return s;
}

string czl_ftos(double f, int digit, int base)
{
    char * s=czl_ftoa(f, digit, base);
    string str;
    if (s!=NULL) {
        str.assign(s);
        delete s;
    }
    return str;
}
// }}}

/**
 * @brief class Message for error/warnning display
 * 
 * This class only contain a ostream pointer for output
 * 
 * @author Zelin Chen
 */
// {{{
//const uint32_t Msg::N=1000;
//uint32_t Msg::nw=0;
//std::vector<std::string> Msg::msg_m;
//std::ofstream Msg::os_m;
//#ifdef PTHREAD
//pthread_mutex_t Msg::lock_m = PTHREAD_MUTEX_INITIALIZER;
//#endif

Msg::Msg()
{
    is_init_m=0;
#ifdef PTHREAD
    pthread_mutex_init(&lock_m, NULL);
#endif
}

Msg::Msg(std::string & file)
{
    file_m = file;
    is_init_m=init(file);
#ifdef PTHREAD
    pthread_mutex_init(&lock_m, NULL);
#endif
}

Msg::Msg(const char file[])
{
    file_m.assign(file);
    is_init_m=init(file);
#ifdef PTHREAD
    pthread_mutex_init(&lock_m, NULL);
#endif
}

Msg::Msg(Msg const & msg)
{
}

Msg::~Msg()
{
//    if (msg_m.size()>0) {
//        for (size_t i=0; i<msg_m.size(); i++) {
//            os_m << msg_m[i] << endl;
//        }
//    }
    os_m.close();
#ifdef PTHREAD
    pthread_mutex_destroy(&lock_m);
#endif
    
}

int Msg::init(const char file[])
{
    if (os_m.is_open()) {
        cerr << "There's a file still open. Do nothing" << std::endl;
        return 0;
    }
    os_m.open(file);
    if ( os_m.fail() ) {
        cerr << "Can't open File " << file << std::endl;
        return 2;
    } else {
        return 0;
    }
}

int Msg::init(string & file)
{
    return init(file.c_str());
}

short Msg::is_init()
{
    return is_init_m;
}

short Msg::is_good()
{
    return good();
}

short Msg::good()
{
    return os_m.good();
}

short Msg::is_fail()
{
    return fail();
}

short Msg::fail()
{
    return os_m.fail();
}

Msg & Msg::print(const char head[], const char msg[], short is_flush, short where)
{
#ifdef PTHREAD
    pthread_mutex_lock(&lock_m);
#endif
    if ( where&0x1) {
        cerr << head << msg << "\n";
        if (is_flush) {
            cerr << std::flush;
        }
    }
    if ( where&0x2) {
        os_m << head << msg << "\n";
        if (is_flush) {
            os_m << std::flush;
        }
    }
#ifdef PTHREAD
    pthread_mutex_unlock(&lock_m);
#endif
    return *this;
}

Msg & Msg::print(string const & head, string const & msg, short is_flush, short where)
{
    return print(head.c_str(), msg.c_str(), is_flush, where);
}

Msg & Msg::print(const char msg[], short is_flush, short where)
{
    return print("", msg, is_flush, where);
}

Msg & Msg::print(const string & msg, short is_flush, short where)
{
    return print("", msg.c_str(), is_flush, where);
}

Msg & Msg::error(const char msg[], const char file[], int line, short is_flush, short where)
{
#ifdef PTHREAD
    pthread_mutex_lock(&lock_m);
#endif
    if ( where&0x1) {
        cerr << "E: " << msg << " (" << file << ":" << line << ")\n"; 
    }
    if ( where&0x2) {
        os_m << "E: " << msg << " (" << file << ":" << line << ")\n"; 
    }
    if (is_flush) {
        cerr << std::flush;
        os_m << std::flush;
    }
#ifdef PTHREAD
    pthread_mutex_unlock(&lock_m);
#endif
    return *this;
}

Msg & Msg::error(stringstream & msg, const char file[], int line, short is_flush, short where)
{
    return error(msg.str().c_str(), file, line, is_flush, where);
}

Msg & Msg::error(stringstream & msg, string & file, int line, short is_flush, short where)
{
    return error(msg.str().c_str(), file.c_str(), line, is_flush, where);
}

Msg & Msg::error(string const & msg, string & file, int line, short is_flush, short where)
{
    return error(msg.c_str(), file.c_str(), line, is_flush, where);
}

Msg & Msg::error(const char msg[], short is_flush, short where)
{
    return print("E: ", msg, is_flush, where);
}

Msg & Msg::error(string const & msg, short is_flush, short where)
{
    return error(msg.c_str(), is_flush, where);
}

Msg & Msg::error(stringstream & msg, short is_flush, short where)
{
    return error(msg.str().c_str(), is_flush, where);
}

Msg & Msg::warn(const char msg[], short is_flush, short where)
{
    return print("W: ", msg, is_flush, where);
}

Msg & Msg::warn(string const & msg, short is_flush, short where)
{
    return warn(msg.c_str(), is_flush, where);
}

Msg & Msg::warn(stringstream & msg, short is_flush, short where)
{
    return warn(msg.str().c_str(), is_flush, where);
}

Msg & Msg::warn(const char msg[], const char file[], int line, short is_flush, short where)
{
#ifdef PTHREAD
    pthread_mutex_lock(&lock_m);
#endif
    if ( where&0x1) {
        cerr << "W: " << msg << "(" << file << ":" << line << ")\n"; 
    }
    if ( where&0x2) {
        os_m << "W: " << msg << "(" << file << ":" << line << ")\n"; 
    }
    if (is_flush) {
        cerr << std::flush;
        os_m << std::flush;
    }
#ifdef PTHREAD
    pthread_mutex_unlock(&lock_m);
#endif
    return *this;
}

Msg & Msg::warn(string const & msg, string & file, int line, short is_flush, short where)
{
    return warn(msg.c_str(), file.c_str(), line, is_flush, where);
}

Msg & Msg::warn(stringstream & msg, string & file, int line, short is_flush, short where)
{
    return warn(msg.str().c_str(), file.c_str(), line, is_flush, where);
}

Msg & Msg::warn(stringstream & msg, const char file[], int line, short is_flush, short where)
{
    return warn(msg.str().c_str(), file, line, is_flush, where);
}

Msg& Msg::operator << (const int i)
{
#ifdef PTHREAD
    pthread_mutex_lock(&lock_m);
#endif
    os_m << i;
//    cerr << i;
#ifdef PTHREAD
    pthread_mutex_unlock(&lock_m);
#endif
    return (*this);
}

Msg& Msg::operator << (const long int i)
{
#ifdef PTHREAD
    pthread_mutex_lock(&lock_m);
#endif
    os_m << i;
//    cerr << i;
#ifdef PTHREAD
    pthread_mutex_unlock(&lock_m);
#endif
    return (*this);
}

Msg& Msg::operator << (const unsigned int i)
{
#ifdef PTHREAD
    pthread_mutex_lock(&lock_m);
#endif
    os_m << i;
//    cerr << i;
#ifdef PTHREAD
    pthread_mutex_unlock(&lock_m);
#endif
    return (*this);
}

Msg& Msg::operator << (const unsigned long int i)
{
#ifdef PTHREAD
    pthread_mutex_lock(&lock_m);
#endif
    os_m << i;
//    cerr << i;
#ifdef PTHREAD
    pthread_mutex_unlock(&lock_m);
#endif
    return (*this);
}

Msg& Msg::operator << (const char i)
{
#ifdef PTHREAD
    pthread_mutex_lock(&lock_m);
#endif
    os_m << i;
//    cerr << i;
#ifdef PTHREAD
    pthread_mutex_unlock(&lock_m);
#endif
    return (*this);
}

Msg& Msg::operator << (const float f)
{
#ifdef PTHREAD
    pthread_mutex_lock(&lock_m);
#endif
    os_m << f;
//    cerr << f;
#ifdef PTHREAD
    pthread_mutex_unlock(&lock_m);
#endif
    return (*this);
}

Msg& Msg::operator << (const double f)
{
#ifdef PTHREAD
    pthread_mutex_lock(&lock_m);
#endif
    os_m << f;
//    cerr << f;
#ifdef PTHREAD
    pthread_mutex_unlock(&lock_m);
#endif
    return (*this);
}

Msg& Msg::operator << (const string & s)
{
#ifdef PTHREAD
    pthread_mutex_lock(&lock_m);
#endif
    os_m << s;
//    cerr << s;
#ifdef PTHREAD
    pthread_mutex_unlock(&lock_m);
#endif
    return (*this);
}

Msg& Msg::operator << (const char *str)
{
#ifdef PTHREAD
    pthread_mutex_lock(&lock_m);
#endif
    os_m << str;
//    cerr << str;
#ifdef PTHREAD
    pthread_mutex_unlock(&lock_m);
#endif
    return (*this);
}

Msg& Msg::operator<< (streambuf* sb )
{
#ifdef PTHREAD
    pthread_mutex_lock(&lock_m);
#endif
    os_m << sb;
//    cerr << sb;
#ifdef PTHREAD
    pthread_mutex_unlock(&lock_m);
#endif
    return (*this);
}

Msg& Msg::operator << (ostream& (*pf)(ostream&))
{
#ifdef PTHREAD
    pthread_mutex_lock(&lock_m);
#endif
    os_m << pf;
//    cerr << pf;
#ifdef PTHREAD
    pthread_mutex_unlock(&lock_m);
#endif
    return (*this);
}

Msg& Msg::operator << (ios& (*pf)(ios&))
{
#ifdef PTHREAD
    pthread_mutex_lock(&lock_m);
#endif
    os_m << pf;
//    cerr << pf;
#ifdef PTHREAD
    pthread_mutex_unlock(&lock_m);
#endif
    return (*this);
}

Msg& Msg::operator << (ios_base& (*pf)(ios_base&))
{
#ifdef PTHREAD
    pthread_mutex_lock(&lock_m);
#endif
    os_m << pf;
//    cerr << pf;
#ifdef PTHREAD
    pthread_mutex_unlock(&lock_m);
#endif
    return (*this);
}

Msg& Msg::endl(Msg & msg)
{
#ifdef PTHREAD
    pthread_mutex_lock(&lock_m);
#endif
    os_m << std::endl;
//    cerr << endl;
#ifdef PTHREAD
    pthread_mutex_unlock(&lock_m);
#endif
    return (*this);
}

Msg& Msg::flush(Msg & msg)
{
#ifdef PTHREAD
    pthread_mutex_lock(&lock_m);
#endif
    os_m << std::flush;
//    cerr << std::flush;
#ifdef PTHREAD
    pthread_mutex_unlock(&lock_m);
#endif
    return (*this);
}

Msg& Msg::flush()
{
#ifdef PTHREAD
    pthread_mutex_lock(&lock_m);
#endif
    os_m.flush();
//    cerr.flush();
#ifdef PTHREAD
    pthread_mutex_unlock(&lock_m);
#endif
    return (*this);
}
// }}}

/* class StringUtility
 */
// {{{
int StringUtility::find(string const & line, char sep, int pos, string &out)
{
    if ( &out == &line) {
        runtime_error e("In Func 'StringUtility::find': input and output are same.");
        throw e;
    }
    int i;
    i = line.find(sep, pos);
    out = line.substr(pos, i-pos);
    if (i!=string::npos) i++;
    return i;
}

int StringUtility::find(string const & line, string const & sep, int pos, string &out)
{
    if ( &out == &line) {
        runtime_error e("In Func 'StringUtility::find': input and output are same.");
        throw e;
    }
    int i;
    int l = sep.size();
    i = line.find(sep, pos);
    out = line.substr(pos, i-pos);
    if (i!=string::npos) i+=l;
    return i;
}

int StringUtility::find_first_of(string const & line, string const & sep, int pos, string &out)
{
    if ( &out == &line) {
        runtime_error e("In Func 'StringUtility::find_first_of': input and output are same.");
        throw e;
    }
    int i;
    i = line.find_first_of(sep, pos);
    out = line.substr(pos, i-pos);
    if (i!=string::npos) i++;
    return i;
}

int StringUtility::skip(string & line, int pos, string sep, int n)
{
    int i = pos;
    while (n>0) {
        i = line.find_first_of(sep, i);
        if (i==string::npos) break;
        i++;
        n--;
    }
    return i;
}

int StringUtility::trim_left(string & str, const char trim)
{
    int i;
    for (i=0; i<str.size(); i++) {
        if (str[i]!=trim) break;
    }
    str = str.erase(0, i);
    return 0;
}

int StringUtility::trim_left(string & str, string const & trim)
{
    int i, k;
    if (trim.size()<10 && str.size()<1024 && trim.size()*str.size()<512) {
        for (i=0; i<str.size(); i++) {
            for (k=0; k<trim.size(); k++) {
                if (str[i]==trim[k]) break;
            }
            if (k>=trim.size()) break;
        }
    } else {
        bool in[256] = {false};
        for (k=0; k<trim.size(); k++) {
            in[ trim[k] ] = true;
        }
        for (i=0; i<str.size(); i++) {
            if (!in[str[i]]) break;
        }
    }
    str = str.erase(0, i);
    return 0;
}

int StringUtility::trim_right(string & str, const char trim)
{
    int i;
    for (i=str.size()-1; i>=0; i--) {
        if (str[i]!=trim) break;
    }
    str = str.erase(i+1);
    return 0;
}

int StringUtility::trim_right(string & str, string const & trim)
{
    int i, k;
    if (trim.size()<10 && str.size()<1024 && trim.size()*str.size()<512) {
        for (i=str.size()-1; i>=0; i--) {
            for (k=0; k<trim.size(); k++) {
                if (str[i]==trim[k]) break;
            }
            if (k>=trim.size()) break;
        }
    } else {
        bool in[256] = {false};
        for (k=0; k<trim.size(); k++) {
            in[ trim[k] ] = true;
        }
        for (i=str.size()-1; i>=0; i--) {
            if (!in[str[i]]) break;
        }
    }
    str = str.erase(i+1);
    return 0;
}

int StringUtility::trim(string & str, const char trim)
{
    int i;
    for (i=0; i<str.size(); i++) {
        if (str[i]!=trim) break;
    }
    int j;
    for (j=str.size()-1; j>i; j--) {
        if (str[j]!=trim) break;
    }
    j++;
    if (j<=i) str.clear();
    else str = str.substr(i, j-i);
    return 0;
}

int StringUtility::trim(string & str, string const & trim)
{
    int i, j, k;
    if (trim.size()<10 && str.size()<1024 && trim.size()*str.size()<512) {
        for (i=0; i<str.size(); i++) {
            for (k=0; k<trim.size(); k++) {
                if (str[i]==trim[k]) break;
            }
            if (k>=trim.size()) break;
        }
        for (j=str.size()-1; j>i; j--) {
            for (k=0; k<trim.size(); k++) {
                if (str[j]==trim[k]) break;
            }
            if (k>=trim.size()) break;
        }
        j++;
    } else {
        bool in[256] = {false};
        for (k=0; k<trim.size(); k++) {
            in[ trim[k] ] = true;
        }
        for (i=0; i<str.size(); i++) {
            if (!in[str[i]]) break;
        }
        for (j=str.size()-1; j>i; j--) {
            if (!in[str[j]]) break;
        }
        j++;
    }
    if (j<=i) str.clear();
    else str = str.substr(i, j-i);
    return 0;
}

int StringUtility::erase_all(string &str, const string & e)
{
    if (e.empty()) return 0;
    int i, j, k;
    j=0;
    if (e.size()<10 && str.size()<1024 && e.size()*str.size()<512) {
        for (i=0; i<str.size(); i++) {
            for (k=0; k<e.size(); k++) {
                if (str[i]==e[k]) break;
            }
            if (k>=e.size()) {
                if (i!=j) str[j] = str[i];
                j++;
            }
        }
    } else {
        bool in[256] = {false};
        for (k=0; k<e.size(); k++) {
            in[ e[k] ] = true;
        }
        for (i=0; i<str.size(); i++) {
            if ( ! in[ str[i] ] ) {
                if (i!=j) str[j] = str[i];
                j++;
            }
        }
    }
    str.resize(j);
    return 0;
}

int StringUtility::split(const string & str, const char sep, vector<string> & out_v)
{
    int i=0;
    string a;
    while (i!=string::npos) {
        i = StringUtility::find(str, sep, i, a);
        out_v.push_back(a);
    }
    return 0;
}

int StringUtility::split(const string & str, const string & sep, vector<string> & out_v)
{
    int i=0;
    string a;
    while (i!=string::npos) {
        i = StringUtility::find(str, sep, i, a);
        out_v.push_back(a);
    }
    return 0;
}


int StringUtility::split_first_of(const string & str, const string & sep, vector<string> & out_v)
{
    int i=0;
    string a;
    while (i!=string::npos) {
        i = StringUtility::find_first_of(str, sep, i, a);
        out_v.push_back(a);
    }
    return 0;
}

int StringUtility::split(const string & str, char sep, char trim, vector<string> & out_v)
{
    int i=0;
    string a;
    while (i!=string::npos) {
        i = StringUtility::find(str, sep, i, a);
        StringUtility::trim(a, trim);
        out_v.push_back(a);
    }
    return 0;
}

int StringUtility::split(const string & str, string const & sep, char trim, vector<string> & out_v)
{
    int i=0;
    string a;
    while (i!=string::npos) {
        i = StringUtility::find(str, sep, i, a);
        StringUtility::trim(a, trim);
        out_v.push_back(a);
    }
    return 0;
}

int StringUtility::split_first_of(const string & str, string const & sep, char trim, vector<string> & out_v)
{
    int i=0;
    string a;
    while (i!=string::npos) {
        i = StringUtility::find_first_of(str, sep, i, a);
        StringUtility::trim(a, trim);
        out_v.push_back(a);
    }
    return 0;
}

void StringUtility::print(const string & str, int indent, int width, ostream & os)
{
    int col = width-indent;
    for (int i=0; i<str.size(); i++) {
        if (i==0) {
            for (int j=0; j<indent; j++) os << ' ';
        } else if (i%col==0) {
            os << "\n";
            for (int j=0; j<indent; j++) os << ' ';
        }
        os << str[i];
    }
}

void StringUtility::to_upper(string & str)
{
    for (int i=0; i<str.size(); i++) {
        str[i] = toupper(str[i]);
    }
}

void StringUtility::to_upper_copy(string const & str, string & out)
{
    out.resize(str.size());
    for (int i=0; i<str.size(); i++) {
        out[i] = toupper(str[i]);
    }
}

void StringUtility::to_lower(string & str)
{
    for (int i=0; i<str.size(); i++) {
        str[i] = tolower(str[i]);
    }
}

void StringUtility::to_lower_copy(string const & str, string & out)
{
    out.resize(str.size());
    for (int i=0; i<str.size(); i++) {
        out[i] = tolower(str[i]);
    }
}

// print
// {{{
string & StringUtility::mar_print(string const &s, int w, int lm, int rm, string & out)
// {{{
{
    if (lm+rm>=w) return out;

    int j=0;
    int l = w-lm-rm;
    for (int i=0; i<s.size(); ) {
        for (j=0; j<lm; j++) out.push_back(' ');

        j=0;
        while ( j<l && i<s.size() ) {
            out.push_back(s[i]);
            i++;
            j++;
        }

        if (j==l) {
            for (j=0; j<rm; j++) out.push_back(' ');
            out.push_back('\n');
        }
    }
    return out;
}
// }}}

ostream & StringUtility::mar_print(string const &s, int w, int lm, int rm, ostream & out)
// {{{
{
    if (lm+rm>=w) return out;

    int j=0;
    int l = w-lm-rm;
    for (int i=0; i<s.size(); ) {
        for (j=0; j<lm; j++) out << ' ';

        j=0;
        while ( j<l && i<s.size() ) {
            out << s[i];
            i++;
            j++;
        }

        if (j==l) {
            for (j=0; j<rm; j++) out << ' ';
            out << '\n';
        }
    }
    
    return out;
}
// }}}
// }}}

// }}}

// class Opt
// {{{
Opt::Opt(char const * names_str, char const * type, char const * desc,
        char const * value, int n, bool is_required)
{
    StringUtility::split(string(names_str), ',', ' ', name_v_m);
    type_m   = type;
    desc_m   = desc;
    StringUtility::trim(value_m, " \t");
    n_m      = n;
    is_required_m = is_required;
    if (!(value==NULL) && strlen(value)>0) {
        value_m  = value;
    }
}

Opt::Opt(string const & names_str, string const & type, 
        string const & desc, string const & value,
        int n, bool is_required)
{
    StringUtility::split(names_str, ',', ' ', name_v_m);
    type_m   = type;
    desc_m   = desc;
    value_m  = value;
    StringUtility::trim(value_m, " \t");
    n_m      = n;
    is_set_m = false;
    is_required_m = is_required;
}

Opt::Opt(vector<string> const & name_v, string const & type,
        string const & desc, string const & value, int n,
        const bool is_required)
{
    name_v_m = name_v;
    type_m   = type;
    desc_m   = desc;
//  p_m    = p;
    value_m  = value;
    StringUtility::trim(value_m, " \t");
    n_m      = n;
    is_set_m = false;
    is_required_m = is_required;
}

void Opt::print_usage(ostream & os, int indent)
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
    stringstream type_ss;
    type_ss << "[" << type_m << "]";
    for (size_t i=1; i<n; i++) {
        os << h << "--" << std::setw(max_name_ol_m) << std::left << name_v_m[i];
        if (n_m==0) {
            os << setw(16) << " " << " : ";
        } else {
            os << setw(16) << std::left << type_ss.str() << " : ";
        }
        os << "\n";
    }
    os << h << "--" << std::setw(max_name_ol_m) << std::left << name_v_m[0];
    if (n_m==0) {
        os << setw(16) << " " << " : ";
    } else {
        os << setw(16) << std::left << type_ss.str() << " : ";
    }
    if (is_required_m) {
        os << "*";
    } else {
        os << " ";
    }
    os << "\n";
    for (i=indent; i<indent+pre_m; i++) h[i] = ' ';
    h[i] = '\0';
    int m0 = 0, m1 = m;
    while (m0<desc_m.size()) {
        os << h << desc_m.substr(m0, m) << "\n";
        m0 += m;
    }
}

void Opt::print_config(ostream & os, int indent)
{
    int i;
    char h[indent+1];
    for (i=0; i<indent; i++) {
        h[i] = ' ';
    }
    h[i] = '\0';
    os << h << std::setw(max_name_ol_m) << std::left << name_v_m[0] << " = " << value_m;
}

// }}}

// class Opts
// {{{
Opts::Opts() 
{
    create<OptS>("config,conf,c", "FILE", "configure file", "config", 1, false);
    create<OptS>("version,ver,v", "", "print the version of this program", "", 0, false);
    create<OptS>("help,h", "", "print the usage", "", 0, false);
}

int Opts::parse(int argc, char** argv)
// {{{
{
    argc_m = argc;
    argv_m = argv;
    int i, id, j, k;
    for (k=1; k<argc; k++) {
        string name;
        string value;
        if ( argv[k][0]=='-' ) {
            if ( argv[k][1]=='-' ) {
                name = argv[k]+2;
            } else {
                name = argv[k]+1;
            }
            StringUtility::trim(name, " \t");
            for (i=0; i<name.size(); i++) {
                if (name[i]=='-') name[i]='_';
            }
            if ( opt_sl_m.find(name) == opt_sl_m.end() ) { /// No this option
                cerr << "E: Option '" << name << "' is not availabe. ";
                cerr << CZL_DBG_INFO << endl;
                return 1;
            } else {
                Opt *p = opt_v_m[ opt_sl_m[name] ];
                if (p->n_m>0) {
                    k++;
                    p->value_m = argv[k];
                    if (p->value_m.empty()) {
                        cerr << "E: Option '" << name 
                            << "' don't have a value.";
                        cerr << CZL_DBG_INFO << endl;
                        return 4;
                    }
                    p->str_to_value();
                }
                if (p->is_set() ) {
                    cerr << "W: Option '" << name << "' occure multipel times.";
                    cerr << CZL_DBG_INFO << endl;
                } else {
                    p->is_set_m = true;
                }
                if ( p->name_v_m[0]=="config") {
                    string conf_file = dynamic_cast<OptS*>(p)->get_value_ref();
                    ifstream fs(conf_file.c_str());
                    if (fs.fail()) {
                        cerr << "E: Can't open configure file '" << conf_file;
                        cerr << "' " << CZL_DBG_INFO << endl;
                        return 2;
                    } else {
                        string str;
                        while (!fs.eof()) {
                            getline(fs, str);
                            StringUtility::trim(str, " \t");

                            if (str.empty()) { continue; }

                            if (str[0]=='[') { continue; }
                            if (str[0]=='#') { continue; }

                            string name1;
                            i=StringUtility::find(str, '=', 0, name1);

                            if (name1.empty()) continue;
                            StringUtility::trim(name1, " \t");
                            if ( opt_sl_m.find(name1)==opt_sl_m.end() ) {
                                cerr << "E: Option '" << name1 << "' is not availabe. ";
                                cerr << CZL_DBG_INFO << endl;
                                return 3;
                            } else {
                                id = opt_sl_m[name1];
                                Opt *p1 = opt_v_m[ id ];
                                p1->value_m=str.substr(i); // set value_m first
                                StringUtility::trim(p1->value_m, " \t");
                                if (p1->is_set() ) {
                                    cerr << "W: Option '" << name1 << "' occure multipel times.";
                                    cerr << CZL_DBG_INFO << endl;
                                } else {
                                    p1->is_set_m = true;
                                }
                                p1->str_to_value();
                            }
                        }
                        fs.close();
                    }
                } else if ( p->name_v_m[0]=="version") {
                    cout << dynamic_cast<OptS*>(p)->get_value_ref() << endl;
                    exit(0);
                } else if ( p->name_v_m[0]=="help") {
                    print_usage(std::cout);
                    exit(0);
                    return 0;
                }
            }
        }
    }
    if ( !is_all_set() ) return 4;
    return 0;
}
// }}}

bool Opts::is_all_set()
{
    int r=0;
    for (vector<Opt*>::iterator it=opt_v_m.begin(); it!=opt_v_m.end(); it++) {
        if ( (*it)->is_required_m && !(*it)->is_set_m) {
            cerr << "Required option '" << (*it)->name_v_m[0] << "' is not set.\n";
            r++;
        }
    }
    cerr << endl;
    return r==0;
}

int Opts::push_back(Opt * p)
{
    try {
        opt_v_m.push_back(p);
    } catch (bad_alloc & ba) {
        std::cerr << "bad_alloc caught: " << ba.what() << CZL_DBG_INFO << '\n';
    }
    for (int i=0; i<p->name_v_m.size(); i++) {
        if ( opt_sl_m.find(p->name_v_m[i])!=opt_sl_m.end() ) {
            cerr << "E: Duplicate option: '" << p->name_v_m[i] << "'" << CZL_DBG_INFO << endl;
            return 1;
        }
        opt_sl_m[ p->name_v_m[i] ] = opt_v_m.size()-1;
    }
    return 0;
}

void Opts::print_usage(ostream & s)
{
    s << "Usage: " << argv_m[0] << " <OPTIONs>\n";
    for (int i=0; i<opt_v_m.size(); i++) {
        opt_v_m[i]->print_usage(s);
        s << "\n";
    }
    s.flush();
}

void Opts::print_config(ostream & s)
{
    for (int i=0; i<opt_v_m.size(); i++) {
        Opt *p = opt_v_m[i];
        if (p->is_set_m) {
            p->print_config(s);
            s << "\n";
        }
    }
    s.flush();
}

int Opts::str_to_value_iv(string const & str, vector<int> & v)
{
    int i=0;
    string str1;
    while (i!=string::npos) {
        StringUtility::find_first_of(str, ",;", i, str1);
        v.push_back(atoi(str1.c_str()));
    }
    return 0;
}

int Opts::str_to_value_lv(string const & str, vector<long> & v)
{
    int i=0;
    string str1;
    while (i!=string::npos) {
        StringUtility::find_first_of(str, ",;", i, str1);
        v.push_back(atol(str1.c_str()));
    }
    return 0;
}

int Opts::str_to_value_fv(string const & str, vector<float> & v)
{
    int i=0;
    string str1;
    while (i!=string::npos) {
        StringUtility::find_first_of(str, ",;", i, str1);
        v.push_back(atof(str1.c_str()));
    }
    return 0;
}

int Opts::str_to_value_lfv(string const & str, vector<double> & v)
{
    int i=0;
    string str1;
    while (i!=string::npos) {
        StringUtility::find_first_of(str, ",;", i, str1);
        v.push_back(atof(str1.c_str()));
    }
    return 0;
}

int Opts::str_to_value_sv(string const & str, vector<string> & v)
{
    int i=0;
    string str1;
    while (i!=string::npos) {
        StringUtility::find_first_of(str, ",;", i, str1);
        StringUtility::trim(str1, " \t");
        v.push_back(str1);
    }
    return 0;
}

Opts::~Opts()
{
    for (vector<Opt*>::iterator it=opt_v_m.begin(); it!=opt_v_m.end(); it++) {
        if (*it != NULL) delete *it;
    }
}
// }}}

// class File
// {{{
bool File::exists_file(string const &fname)
{
    ifstream fin(fname.c_str());
    if (fin.fail()) return false;
    struct stat st;
    if ( stat(fname.c_str(), &st) ) { throw; }
    if ( S_ISREG(st.st_mode) ) return true;
    else return false;
}

bool File::exists_dir(string const &fname)
{
    struct stat st;
    if ( stat(fname.c_str(), &st) ) { return false; }
    if ( S_ISDIR(st.st_mode) ) return true;
    else return false;
}

int File::create_dir(string const &fname)
{
    return mkdir(fname.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
}
// }}}


bool is_be()
{
    uint8_t a[4] = {0x01, 0x0, 0x0, 0x0};
    uint32_t b;
    memcpy(&b, a, sizeof(uint32_t));
    return b!=1;
}

int machine_bit()
{
    return sizeof(void*)*8;
}

string str_time(time_t const & t)
{
    char s[128];
    strftime(s, 128, "%c", localtime(&t));
    string str(s);
    return str;
}

string str_difftime(time_t const & end, time_t const & begin)
{
    int i;
    double d = difftime(end, begin);
    int m1 = 60;
    int h1 = m1*60;
    int d1 = h1*24;
    int day = static_cast<long>(d)/d1;
    i = static_cast<long>(d)%d1;
    int hour = i/h1;
    i = i%h1;
    int minute = i/m1;
    int second = i%m1;
    stringstream ss;
    ss << day << "d " << hour << "h " << minute << "m " << second << "s"; 
    return ss.str();
}

};
