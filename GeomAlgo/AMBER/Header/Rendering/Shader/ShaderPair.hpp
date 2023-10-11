#ifndef __ENGINE_SHADER_PAIR__
#define __ENGINE_SHADER_PAIR__

#include <iostream>
namespace Ge
{
	class ShaderPair
	{
	public:
		ShaderPair(const std::string & f, const std::string & v,bool b,bool m = true,bool t = false,int cm = 1);
		ShaderPair();
		std::string Frag;
		std::string Vert;
		bool back;
		bool multiSampling = true;
		bool transparency = false;
		int cullMode = 1;
	};
}

#endif //__ENGINE_SHADER_PAIR__