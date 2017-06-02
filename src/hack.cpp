/*	Hacks :
		unusual codes that may help in contests.
*/

//	long long formatting hack :

#ifdef WIN32
	#define LL "%I64d"
#else
	#define LL "%lld"
#endif

//	Optimizing hack :

#pragma GCC optimize ("O3")
#pragma GCC optimize ("whole-program")

//	Stack hack :

//	C++

#pragma comment(linker, "/STACK:36777216")

//	G++

int __size__ = 256 << 20;	//	256MB
char *__p__ = (char*)malloc(__size__) + __size__;
__asm__("movl %0, %%esp\n" :: "r"(__p__));
