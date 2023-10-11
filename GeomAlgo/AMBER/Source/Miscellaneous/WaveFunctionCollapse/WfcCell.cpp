#include "WfcCell.hpp"

WfcCell::WfcCell()
{

}

WfcCell::WfcCell(std::vector<int> & all_posibility)
{
    std::copy(all_posibility.begin(), all_posibility.end(), m_possibilites.begin());
}

void WfcCell::copy(std::vector<int>& all_posibility)
{
    std::copy(all_posibility.begin(), all_posibility.end(), m_possibilites.begin());
}

void WfcCell::collapse()
{
    if (this->m_possibilites.size() >= 1)
    {
        int index = std::rand() % this->m_possibilites.size();
        int t = this->m_possibilites.at(index);
        this->m_possibilites.clear();
        this->m_possibilites.push_back(t);
    }
    else
    {
        this->m_possibilites.clear();
        this->m_possibilites.push_back(-1);
    }

    this->m_value = this->m_possibilites.at(0);
    this->m_collapsed = true;

    return;
}

int WfcCell::get_entropy()
{
    return this->m_possibilites.size();
}

int WfcCell::get_value()
{
    return m_value;
}

void WfcCell::update(int value, Direction d)
{

}

bool WfcCell::is_collapsed()
{
	return this->m_collapsed;
}
