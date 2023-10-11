#ifndef __SHADER_UTIL__
#define __SHADER_UTIL__

#include <vector>

class ShaderUtil
{
public:   
    static void CalcWorkSize(int length, int * x, int * y, int * z)
    {
        int GROUP_SIZE = 256;
        int MAX_DIM_GROUPS = 256;
        int MAX_DIM_THREADS = (GROUP_SIZE * MAX_DIM_GROUPS);
        int MAX_DIM_THREADS_THREADS = (MAX_DIM_THREADS * MAX_DIM_GROUPS);
        if (length <= MAX_DIM_THREADS)
        {
            *x = (length - 1) / GROUP_SIZE + 1;
            *y = *z = 1;
        }
        else if (length <= MAX_DIM_THREADS_THREADS)
        {
            *x = MAX_DIM_GROUPS;
            *y = (length - 1) / MAX_DIM_THREADS + 1;
            *z = 1;
        }
        else
        {
            *x = *y = MAX_DIM_GROUPS;
            *z = (length - 1) / MAX_DIM_THREADS_THREADS + 1;
        }
    }

    static std::vector<int> DecomposeFirstFactors(int n)
    {
        std::vector<int> facteursPremiers = std::vector<int>();

        while (n % 2 == 0)
        {
            facteursPremiers.push_back(2);
            n /= 2;
        }

        for (int i = 3; i <= sqrt(n); i += 2)
        {
            while (n % i == 0)
            {
                facteursPremiers.push_back(i);
                n /= i;
            }
        }

        if (n > 2)
        {
            facteursPremiers.push_back(n);
        }

        return facteursPremiers;
    }
};


#endif // !__SHADER_UTIL__
