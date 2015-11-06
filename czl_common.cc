/* 
 * czl_common.h
 * @brief Common Functions
 * @author CZL
 */
#include "czl_common.h"

using namespace std;

namespace czl_bio {

/*
 * string -- number convertor
 */
// {{{
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
    n = f;

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

Msg::Msg(Msg & msg)
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

void Msg::error(const char msg[], const char file[], int line, short is_flush)
{
#ifdef PTHREAD
    pthread_mutex_lock(&lock_m);
#endif
    cerr << "E: " << msg << " (" << file << ":" << line << ")\n"; 
    os_m << "E: " << msg << " (" << file << ":" << line << ")\n"; 
    if (is_flush) {
        cerr << std::flush;
        os_m << std::flush;
    }
    os_m.close();
#ifdef PTHREAD
    pthread_mutex_unlock(&lock_m);
#endif
    exit(1);
}

void Msg::error(stringstream & msg, const char file[], int line, short is_flush)
{
    error(msg.str().c_str(), file, line, is_flush);
}

void Msg::error(stringstream & msg, string & file, int line, short is_flush)
{
    error(msg.str().c_str(), file.c_str(), line, is_flush);
}

void Msg::error(string & msg, string & file, int line, short is_flush)
{
    error(msg.c_str(), file.c_str(), line, is_flush);
}

void Msg::error(const char msg[], short is_flush)
{
#ifdef PTHREAD
    pthread_mutex_lock(&lock_m);
#endif
    
//    for (size_t i=0; i<msg_m.size(); i++) {
//        os_m << msg_m[i] << endl;
//    }
//    msg_m.clear();
    cerr << "E: " << msg << "\n"; 
    os_m << "E: " << msg << "\n"; 
    if (is_flush) {
        cerr << std::flush;
        os_m << std::flush;
    }
    os_m.close();
#ifdef PTHREAD
    pthread_mutex_unlock(&lock_m);
#endif
    exit(1);
}

void Msg::error(string & msg, short is_flush)
{
    error(msg.c_str(), is_flush);
}

void Msg::error(stringstream & msg, short is_flush)
{
    error(msg.str().c_str(), is_flush);
}

void Msg::warn(const char msg[], short is_flush)
{
//    string msg1 = "W ";
//    char a[128];
//    sprintf(a, "%i", nw);
//    nw++;
//    msg1 = msg1+ a + ": " + msg;
#ifdef PTHREAD
    pthread_mutex_lock(&lock_m);
#endif
//    msg_m.push_back(msg1);
//    if (msg_m.size()>=N) {
//        for (size_t i=0; i<msg_m.size(); i++) {
//            os_m << msg_m[i] << endl;
//        }
//        msg_m.clear();
//    }
//    if (nw%N == 0) {
//        std::cerr << "More than " << nw << " Warnnings." << endl;
//    //  std::cout << "More than " << nw << " Warnnings. Continue? [Y|n]: ";
//        char c;
//    //  std::cin >> c;
//    //  if (c=='N' || c=='n' || c=='\n') {
//    //      exit(1);
//    //  }
//    }
    cerr << "W: " << msg << "\n"; 
    os_m << "W: " << msg << "\n"; 
    if (is_flush) {
        cerr << std::flush;
        os_m << std::flush;
    }
#ifdef PTHREAD
    pthread_mutex_unlock(&lock_m);
#endif
}

void Msg::warn(string & msg, short is_flush)
{
    warn(msg.c_str(), is_flush);
}

void Msg::warn(stringstream & msg, short is_flush)
{
    warn(msg.str().c_str(), is_flush);
}

void Msg::warn(const char msg[], const char file[], int line, short is_flush)
{
#ifdef PTHREAD
    pthread_mutex_lock(&lock_m);
#endif
    cerr << "W: " << msg << "(" << file << ":" << line << ")\n"; 
    os_m << "W: " << msg << "(" << file << ":" << line << ")\n"; 
    if (is_flush) {
        cerr << std::flush;
        os_m << std::flush;
    }
#ifdef PTHREAD
    pthread_mutex_unlock(&lock_m);
#endif
}

void Msg::warn(string & msg, string & file, int line, short is_flush)
{
    warn(msg.c_str(), file.c_str(), line, is_flush);
}

void Msg::warn(stringstream & msg, string & file, int line, short is_flush)
{
    warn(msg.str().c_str(), file.c_str(), line, is_flush);
}

void Msg::warn(stringstream & msg, const char file[], int line, short is_flush)
{
    warn(msg.str().c_str(), file, line, is_flush);
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

    /* class Token
     */
    // {{{
    int StringUtility::get_next(string & line, char sep, int pos, string sep, string &out)
    {
        if ( &out == &line) {
            runtime_error e("In Func 'StringUtility::get_next': input and output are same.")
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

    int StringUtility::trim_right(string & str, const char trim)
    {
        int i;
        for (i=str.size()-1; i>=0; i--) {
            if (str[i]!=trim) break;
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
        for (j=0; j<str.size(); j++) {
            if (str[j]!=trim) break;
        }
        if (j<=i) str.clear();
        else str = str.substr(i, j-i);
        return 0;
    }

    int StringUtility::str_to_vector(const string & str, const char sep, vector<string> & out_v)
    {
        int i=0;
        string a;
        while (i!=string::npos) {
            i = this->get_next(str, sep, i, a);
            out_v.push_back(a);
        }
    }

    int StringUtility::str_to_vector(const string & str, const char sep, const char trim, vector<string> & out_v);
    {
        int i=0;
        string a;
        while (i!=string::npos) {
            i = this->get_next(str, sep, i, a);
            trim(a, trim);
            out_v.push_back(a);
        }
    }
    // }}}
};
