#include <iostream>
#include <utility>

struct Node
{
    int data;
    Node* next;
    Node* prev;
    Node(int value) : data(std::move(value)),
        next(nullptr),
        prev(nullptr)
    {}
};

class circular_doubly_linked_list
{
    
    

public:
    Node* head, * tail;

    circular_doubly_linked_list() : head(nullptr),
        tail(nullptr)
    {}
    //copy constructor
    circular_doubly_linked_list(const circular_doubly_linked_list&);

    circular_doubly_linked_list(circular_doubly_linked_list&&) noexcept;

    ~circular_doubly_linked_list();

    void append_node(int value);
    void delete_node(int value);
    void delete_node(Node* value);

    friend void swap(circular_doubly_linked_list& lhs, circular_doubly_linked_list& rhs)
    {
        std::swap(lhs.head, rhs.head);
    }


private:

    struct Node* search(int value)
    {
        Node* node = head;
        while (node->next != head)
        {
            if (node->data == value)
            {
                return node;
            }
            node = node->next;
        }
        if (node->data == value)
        {
            return node;
        }
        std::cerr << "No such element in the list \n";
        return nullptr;
    }

    void print_list(std::ostream& os = std::cout) const
    {
        Node* tmp = head;
        while (tmp->next != head)
        {
            std::cout << tmp->data << ' ';
            tmp = tmp->next;
        }
        std::cout << tmp->data << '\n';
    }
};

circular_doubly_linked_list::circular_doubly_linked_list(const circular_doubly_linked_list& cdll)
{
    if (cdll.head == nullptr)
    {
        head = tail = nullptr;
    }
    else
    {
        head = new Node(cdll.head->data);
        Node* curr = head;
        Node* tmp = head;
        Node* obj_curr = cdll.head;

        while (obj_curr->next != cdll.head)
        {
            curr->next = new Node(obj_curr->next->data);
            obj_curr = obj_curr->next;
            curr = curr->next;
            curr->prev = tmp;
            tmp = tmp->next;
        }
        tail = curr;
        curr->next = head;
        head->prev = curr;
    }
}

circular_doubly_linked_list::circular_doubly_linked_list(circular_doubly_linked_list&& cdll) noexcept
{
    head = tail = nullptr;
    swap(*this, cdll);
}

void circular_doubly_linked_list::append_node(int value)
{
    Node* node = new Node(std::move(value));

    if (head == nullptr)
    {
        node->next = node;
        node->prev = node;
        head = node;
        tail = node;
    }

    tail = head->prev;
    tail->next = node;
    node->prev = tail;
    node->next = head;
    head->prev = node;
    tail = node;
}

void circular_doubly_linked_list::delete_node(int value)
{
    Node* node = search(value);
    if (node == nullptr)
    {
        std::cerr << "No such value in the list\n";
        return;
    }
    else
    {
        Node* tmp = head;
        Node* tail = head->prev;
        if (tmp == node)
        {
            tail->next = tmp->next;
            tmp->prev->next->prev = tail;
            head = tail->prev;
            delete tmp;
            return;
        }
        else if (tail == node)
        {
            Node* curr = tail;
            tmp = tail->prev;
            tmp->next = head;
            head->prev = tmp;
            tail = tmp;
            delete curr;
            return;
        }
        else
        {
            while (tmp->next != head)
            {
                if (tmp == node)
                {
                    tmp->prev->next = tmp->next;
                    tmp->prev->next->prev = tmp->prev;
                    delete tmp;
                    return;
                }
                tmp = tmp->next;
            }
        }
    }
}
void circular_doubly_linked_list::delete_node(Node* value)
{
    Node* node = value;
    if (node == nullptr)
    {
        std::cerr << "No such value in the list\n";
        return;
    }
    else
    {
        Node* tmp = head;
        Node* tail = head->prev;
        if (tmp == node)
        {
            tail->next = tmp->next;
            tmp->prev->next->prev = tail;
            head = tail->prev;
            delete tmp;
            return;
        }
        else if (tail == node)
        {
            Node* curr = tail;
            tmp = tail->prev;
            tmp->next = head;
            head->prev = tmp;
            tail = tmp;
            delete curr;
            return;
        }
        else
        {
            while (tmp->next != head)
            {
                if (tmp == node)
                {
                    tmp->prev->next = tmp->next;
                    tmp->prev->next->prev = tmp->prev;
                    delete tmp;
                    return;
                }
                tmp = tmp->next;
            }
        }
    }
}

circular_doubly_linked_list::~circular_doubly_linked_list()
{
    if (head)
    {
        Node* tmp = head;
        while (tmp->next != head)
        {
            Node* t = tmp;
            tmp = tmp->next;
            delete t;
        }
        delete tmp;
        head = nullptr;
    }
}