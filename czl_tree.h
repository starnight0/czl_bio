#include "czl_common.h"

namespace czl_bio {
	template <typename T> class Tree {
		class TreeNode {
			TreeNode* get_parent()
			{
				return link_m[0];
			}
			void set_parent(TreeNode* p)
			{
				if (link_m.size()==0) link_m.resize(1);
				link_m[0] = p;
			}
			T & get_data_ref()
			{
				return data_m;
			}
			void set_data(const T & data)
			{
				data_m = data;
			}

			T data_m;
			vector<TreeNode*> link_m; // 0: link to parent, other: link to child
		};
	};
};
