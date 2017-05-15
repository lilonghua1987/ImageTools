#pragma once

#include <exception>

namespace alg
{
	template<class T>
	/*The default sort rule is ascend.*/
	void quickSort(T data[], size_t length, bool asc = true)
	{
		if (data == nullptr) throw(std::exception("invalid input"));
		if (length <= 1) return;

		T k = *data;
		int i = 0, j = length - 1;

		while (i < j)
		{
			if (asc)
			{
				//right -> left
				while (i < j && data[j] >= k) j--;

				if (i < j && data[j] < k)
				{
					std::swap(data[i], data[j]);
					i++;
				}

				//left -> right
				while (i < j && data[i] <= k) i++;

				if (i < j && data[i] > k)
				{
					std::swap(data[i], data[j]);
					j--;
				}
			}
			else
			{
				//right -> left
				while (i < j && data[j] <= k) j--;

				if (i < j && data[j] > k)
				{
					std::swap(data[i], data[j]);
					i++;
				}

				//left -> right
				while (i < j && data[i] >= k) i++;

				if (i < j && data[i] < k)
				{
					std::swap(data[i], data[j]);
					j--;
				}
			}
		}

		if (i > 0) quickSort<T>(data,i + 1,asc);
		if (j < (int)length - 1) quickSort<T>(&data[j + 1], length - j - 1,asc);
	}


	template<class T>
	void shellSort(T data[], size_t length, bool asc = true)
	{ 
		if (data == nullptr) throw(std::exception("Invalid input: the sort data must length than 1 !"));
		if (length < 2) throw(std::exception("Invalid input: the lenght must greater than 1 !"));

		size_t d = length / 2;
		while (d >= 1)
		{
			for (size_t i = d; i < length; i++)
			{
				int temp = data[i];
				int j = (int)(i - d);
				if (asc)
				{
					for (; j >= 0 && (data[j] > temp); j -= d)
					{
						data[j + d] = data[j];
					}
				}
				else
				{
					for (; j >= 0 && (data[j] < temp); j -= d)
					{
						data[j + d] = data[j];
					}
				}
				
				data[j + d] = temp;
			}

			d /= 2;
		}
	}

	template<class T>
	class Node
	{
		public:
			Node()
				: parent(nullptr)
				, left(nullptr)
				, right(nullptr)
				, value(0)
				, visited(false)
			{};

			Node(Node* parent, Node* left, Node* right, T value, bool visited = false)
				:parent(parent)
				, left(left)
				, right(right)
				, value(value)
				, visited(visited)
			{};

			~Node()
			{
				if (parent) delete parent;
				if (left) delete left;
				if (right) delete right;
			};

			Node* find(Node* node)
			{
				if (node == nullptr) throw(std::exception("Invalid input: the node must not be null !"));

				node->visited = true;

				if (node->left == nullptr || node->right == nullptr) return node;

				if (node->parent == nullptr) 
					find(node->left);
				else 
					if (node->parent->right->visited) 
						find(node->left);
					else
                        find(node->right);
			}

			Node* parent;
			Node* left;
			Node* right;
			T value;
			bool visited;

	   private:
			Node& operator = (const Node& node)
			{
				return *this;
			}

	};

	template<class T>
	void heapSort(T data[], size_t length)
	{
		Node<T>* head = nullptr;
		Node<T>* temp = nullptr;

		for (size_t i = 1; i < length; i++)
		{
			if (head == nullptr)
			{
				head = new Node<T>(nullptr, nullptr, nullptr, data[i]);
				temp = head;
			}
			else
			{
				Node<T>* node = new Node<T>(temp, nullptr, nullptr, data[i]);

				if (node && !node->parent->left) node->parent->left = node;
				else if (node) node->parent->right = node;

				if (node && node->parent->left && node->parent->right)
				{
					//if (node->)
				}
			}
		}
	}
}