#include "ShaderPair.hpp"

namespace Ge
{
	ShaderPair::ShaderPair(const std::string & f, const std::string &v, bool b,bool m,bool t,int cm)
	{
		this->Frag = f;
		this->Vert = v;
		this->back = b;
		this->multiSampling = m;
		this->transparency = t;
		this->cullMode = cm;
	}

	ShaderPair::ShaderPair(){}
}