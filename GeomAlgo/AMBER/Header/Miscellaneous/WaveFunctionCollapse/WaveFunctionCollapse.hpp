#ifndef __WAVE_FUNCTION_COLLAPSE__
#define __WAVE_FUNCTION_COLLAPSE__

#include "WfcCell.hpp"
#include <vector>
#include <iostream>
#include <map>

template<int WIDTH, int HEIGHT>
class WaveFunctionCollapse
{
public:
    WaveFunctionCollapse(int** paterne,int sizeX,int sizeY);
    void step();
    void generate();
    bool is_solved();
private:
    WfcCell m_arr[WIDTH * HEIGHT];
    int** m_paterne;
    int m_sizeX;
    int m_sizeY;
};

template<int WIDTH, int HEIGHT>
WaveFunctionCollapse<WIDTH, HEIGHT>::WaveFunctionCollapse(int** paterne, int sizeX, int sizeY)
{
    m_paterne = paterne;
    m_sizeX = sizeX;
    m_sizeY = sizeY;
    std::vector<int> m_allPosibility;
    std::map<int, std::map<Direction, std::vector<int>>> map_possibility;
    for (int x = 0; x < m_sizeX; x++) 
    {
        for (int y = 0; y < m_sizeY; y++) 
        {
            int valeur = m_paterne[x][y];
            if (map_possibility.find(valeur) == map_possibility.end()) 
            {
                map_possibility[valeur] = {};
            }

            // Coordonnées des voisins
            int voisinsX[4] = { x, x, x - 1, x + 1 };
            int voisinsY[4] = { y - 1, y + 1, y, y };

            for (int i = 0; i < 4; i++) 
            {
                if (voisinsX[i] >= 0 && voisinsX[i] < m_sizeX && voisinsY[i] >= 0 && voisinsY[i] < m_sizeY) 
                {
                    Direction direction = (i == 0) ? Direction::D_UP :
                        (i == 1) ? Direction::D_DOWN :
                        (i == 2) ? Direction::D_LEFT :
                        Direction::D_RIGHT;

                    auto& vecteur = map_possibility[valeur][direction];
                    int voisin = m_paterne[voisinsX[i]][voisinsY[i]];

                    if (std::find(vecteur.begin(), vecteur.end(), voisin) == vecteur.end()) 
                    {
                        vecteur.push_back(voisin);
                    }
                }
            }
        }
    }

    for (const auto& entry : map_possibility)
    {
        m_allPosibility.push_back(entry.first);        
    }

    for (int i = 0; i < WIDTH * HEIGHT; i++)
    {
        this->m_arr[i] = WfcCell(m_allPosibility);

    }

    /*for (const auto& entry : map_possibility)
    {
        int key = entry.first;
        std::cout << "Key: " << key << std::endl;
    
        for (const auto& directionEntry : entry.second) 
        {
            Direction direction = directionEntry.first;
            std::cout << "  Direction: " << direction << std::endl;
        
            const std::vector<int>& values = directionEntry.second;
            for (const auto& value : values) 
            {
                std::cout << "    Value: " << value << std::endl;
            }
        }
    }*/
}

template<int WIDTH, int HEIGHT>
void WaveFunctionCollapse<WIDTH, HEIGHT>::step()
{
    std::vector<int> lowest_entropy;

    int current_lowest = INT32_MAX;

    for (int i = 0; i < WIDTH * HEIGHT; i++)
    {
        int cell_entropy = this->m_arr[i].get_entropy();

        if (this->m_arr[i].is_collapsed())
        {
            continue;
        }
        else if (cell_entropy == current_lowest)
        {
            lowest_entropy.push_back(i);
        }
        else if (cell_entropy < current_lowest)
        {
            lowest_entropy.clear();
            lowest_entropy.push_back(i);
            current_lowest = cell_entropy;
        }
    }

    if (lowest_entropy.size() > 0)
    {
        int index = std::rand() % lowest_entropy.size();

        WfcCell& c = this->m_arr[lowest_entropy.at(index)];

        c.collapse();

        int tile = this->m_arr[lowest_entropy.at(index)].get_value();

        // top
        if (lowest_entropy.at(index) - WIDTH >= 0)
        {
            this->m_arr[lowest_entropy.at(index) - WIDTH].update(tile, Direction::D_DOWN);
        }

        // bottom
        if (lowest_entropy.at(index) + WIDTH < WIDTH * HEIGHT)
        {
            this->m_arr[lowest_entropy.at(index) + WIDTH].update(tile, Direction::D_UP);
        }

        // left
        if (lowest_entropy.at(index) - 1 >= 0)
        {
            this->m_arr[lowest_entropy.at(index) - 1].update(tile, Direction::D_RIGHT);
        }

        // right
        if (lowest_entropy.at(index) + 1 < WIDTH * HEIGHT)
        {
            this->m_arr[lowest_entropy.at(index) + 1].update(tile, Direction::D_LEFT);
        }
    }
}

template<int WIDTH, int HEIGHT>
bool WaveFunctionCollapse<WIDTH, HEIGHT>::is_solved()
{
    bool solved = true;

    for (int i = 0; i < WIDTH * HEIGHT; i++)
    {
        if (!this->m_arr[i].is_collapsed())
        {
            solved = false;
            break;
        }
    }

    return solved;
}

template<int WIDTH, int HEIGHT>
void WaveFunctionCollapse<WIDTH, HEIGHT>::generate()
{
    while (!this->is_solved())
    {
        this->step();
        int collapsed_count = 0;

        for (WfcCell& c : this->m_arr)
        {
            if (c.is_collapsed())
            {
                collapsed_count++;
            }
        }

    }

    return;
}

#endif //!__WAVE_FUNCTION_COLLAPSE__