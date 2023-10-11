#ifndef __WFC_CELL__
#define __WFC_CELL__

#include <vector>

enum Direction
{
    D_UP,
    D_DOWN,
    D_LEFT,
    D_RIGHT
};


class WfcCell
{    
public:
    WfcCell();
    WfcCell(std::vector<int>& all_posibility);

    void copy(std::vector<int>& all_posibility);
    void collapse();
    int get_entropy();
    bool is_collapsed();
    int get_value();
    void update(int value, Direction d);

private:
    std::vector<int> m_possibilites;
    bool m_collapsed;
    int m_value;
};

#endif//!__WFC_CELL__