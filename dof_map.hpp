#ifndef MARS_DOF_MAP_HPP
#define MARS_DOF_MAP_HPP

#include <algorithm>

namespace mars {

	class Map {
	public:
		Map(const Integer partition_id,
			const Integer n_partitions)
		: partition_id_(partition_id), n_partitions_(n_partitions)
		{}

		inline Integer n_partitions() const
		{
			return n_partitions_;
		}

		// inline Integer size(const Integer partition_id) const
		// {
		// 	return n_x_partition_[partition_id + 1] - n_x_partition_[partition_id];
		// }

		inline Integer max_id() const
		{
			if(global_id_.empty()) return INVALID_INDEX;

			return *std::max_element(global_id_.begin(), global_id_.end());
		}

		inline void push_global(
			const Integer global_id,
			const Integer owner = INVALID_INDEX)
		{
			Integer local_id = global_id_.size();
			global_id_.push_back(global_id);
			global_to_local_[global_id] = local_id;
			owner_id_.push_back(owner);
		}

		inline void add(
			const Integer local_id,
			const Integer global_id,
			const Integer owner = INVALID_INDEX)
		{
			if(local_id >= global_id_.size()) {
				global_id_.resize(local_id + 1, INVALID_INDEX);
			}

			set_global(local_id, global_id);
			set_owner(local_id, owner);
		}

		inline void set_owner(
			const Integer local_id,
			const Integer owner)
		{
			if(local_id >= owner_id_.size()) {
				owner_id_.resize(local_id + 1, INVALID_INDEX);
			}

			if(owner != INVALID_INDEX) { 
				owner_id_[local_id] = owner;
			}
		}

		inline Integer owner(const Integer local_id) const
		{
			if(local_id >= owner_id_.size() || local_id < 0) return INVALID_INDEX;
			// assert(local_id < owner_id_.size());
			return owner_id_[local_id];
		}

		inline void set_global(const Integer local_id, const Integer global_id)
		{
			assert(local_id < global_id_.size());
			global_id_[local_id] = global_id;

			if(global_id != INVALID_INDEX) {
				global_to_local_[global_id] = local_id;
			}
		}

		inline void insert_local_to_global(const Integer local_id, const Integer global_id)
		{
			if(local_id >= global_id_.size()) {
				global_id_.resize(local_id + 1, INVALID_INDEX);
			}

			global_id_[local_id] = global_id;

			if(global_id != INVALID_INDEX) {
				global_to_local_[global_id] = local_id;
			}
		}

		inline Integer local(const Integer global_id) const
		{
			auto it = global_to_local_.find(global_id);
			if(it == global_to_local_.end()) return INVALID_INDEX;

			return it->second;
		}

		inline Integer global(const Integer local_id) const
		{
			// assert(local_id < global_id_.size());
			if(local_id >= global_id_.size() || local_id < 0) return INVALID_INDEX;

			return global_id_[local_id];
		}


		inline Integer insert_global(const Integer local_id)
		{
			// assert(local_id < global_id_.size());
			if(local_id >= global_id_.size()) {
				global_id_.resize(local_id + 1, INVALID_INDEX);
			}

			return global_id_[local_id];
		}

		inline Integer insert_owner(const Integer local_id)
		{
			// assert(local_id < global_id_.size());
			if(local_id >= owner_id_.size()) {
				owner_id_.resize(local_id + 1, INVALID_INDEX);
			}

			return owner_id_[local_id];
		}

		inline Integer n_owned() const
		{
			Integer ret = 0;
			for(auto n : owner_id_) {
				ret += (n == partition_id());
			}

			return ret;
		}

		inline Integer n_unassigned() const
		{
			Integer ret = 0;
			for(auto n : owner_id_) {
				ret += (n == INVALID_INDEX);
			}

			return ret;
		}

		inline Integer n_without_global() const
		{
			Integer ret = 0;
			for(auto n : global_id_) {
				ret += (n == INVALID_INDEX);
			}

			return ret;
		}

		inline Integer n_owned_without_global() const
		{
			Integer ret = 0;
			for(Integer i = 0; i < global_id_.size(); ++i) {
				if(owner(i) == partition_id()) {
					ret += global(i) == INVALID_INDEX;
				}
			}

			return ret;
		}

		inline Integer n_owned_with_global() const
		{
			Integer ret = 0;
			for(Integer i = 0; i < global_id_.size(); ++i) {
				if(owner(i) == partition_id()) {
					ret += global(i) != INVALID_INDEX;
				}
			}

			return ret;
		}

		inline Integer partition_id() const
		{
			return partition_id_;
		}

		void assign_global_ids(
			const Integer begin,
			const Integer end)
		{
			if(begin == end) return;

			Integer current_id = begin;
			
			Integer count = 0;
			for(Integer i = 0; i < global_id_.size(); ++i) {
				if(owner(i) != partition_id() || global(i) != INVALID_INDEX) continue;
					set_global(i, current_id++);
					count++;
			}

			assert(count == end - begin);
		}

		inline const std::vector<Integer> &global_id() const
		{
			return global_id_;
		}

		void pack_for_global(
			std::vector<Integer> &buffer,
			const bool clear_buffer) const
		{
			if(clear_buffer) {
				buffer.clear();
			}

			buffer.resize(n_partitions() + 1, 0);

			Integer n = n_owned_without_global();
			buffer[partition_id() + 1] = n;

			//Hack
			if(clear_buffer) {
				buffer[0] = max_id() + 1;
			} else {
				buffer[0] = std::max(max_id() + 1, buffer[0]);
			}
		}

		void unpack_for_global(const std::vector<Integer> &buffer)
		{
			Integer range_begin = 0;
			Integer range_end   = 0;

			for(Integer i = 0; i <= partition_id(); ++i) {
				range_begin += buffer[i];
			}

			range_end = range_begin + buffer[partition_id() + 1];

			assign_global_ids(range_begin, range_end);
		}

		void describe(std::ostream &os) const
		{
			os << "--------------------\n";
			os << "partition: " << partition_id_ << "\n";

			for(Integer i = 0; i < global_id_.size(); ++i) {
				os << "(" << i << ", " << global(i) << ") " << owner(i) << std::endl;
			}

			os << "--------------------\n";
		}

		void resize(const Integer n, const Integer default_owner)
		{
			global_id_.resize(n, INVALID_INDEX);
			owner_id_.resize(n, default_owner);
			partitions_.resize(n);
		}

		void add_partition(const Integer local_id, const Integer partition_id)
		{
			if(partitions_.size() <= local_id) {
				partitions_.resize(local_id+1);
			}

			partitions_[local_id].push_back(partition_id);
		}

		template<class Iter>
		void set_partitions(const Integer local_id, const Iter begin, const Iter end)
		{
			if(partitions_.size() <= local_id) {
				partitions_.resize(local_id+1);
			}
			
			assert(std::distance(begin, end) > 0);

			partitions_[local_id].clear();
			partitions_[local_id].insert(partitions_[local_id].end(), begin, end);
			assert(is_unique(partitions_[local_id]));
		}

		bool is_unique(const std::vector<Integer> &vec) const
		{
			for(std::size_t i = 1; i < vec.size(); ++i) {
				if(vec[i] == vec[i-1]) return false;
			}

			return true;
		}

		bool has_partitions(const Integer local_id) const
		{
			if(local_id < partitions_.size()) {
				return !partitions_[local_id].empty();
			}

			return false;
		}

		const std::vector<Integer> &partitions(const Integer local_id) const
		{
			assert(local_id < partitions_.size());
			return partitions_[local_id];
		}

		template<class Iter>
		void intersect_partitions(
			const Iter &begin,
			const Iter &end,
			std::vector<Integer> &out) const
		{
			out.clear();
			if(begin == end) return;
			
			auto it = begin;
			auto dist = std::distance(begin, end);
			if(dist == 1) {
				out = partitions_[*it];
				assert(is_unique(out));
				return;
			}

			auto first_it  = *it++;
			auto second_it = *it++;

			std::set_intersection(
				partitions_[first_it].begin(),  partitions_[first_it].end(),
				partitions_[second_it].begin(), partitions_[second_it].end(),
				std::back_inserter(out));

			std::vector<Integer> temp = out;
			for(; it != end; ++it) {
				out.clear();
				std::set_intersection(
					temp.begin(), temp.end(),
					partitions_[*it].begin(), partitions_[*it].end(),
					std::back_inserter(out)
				);

				temp = out;
			}

			assert(is_unique(out));
		}

		bool is_valid(const bool check_partitions = true) const
		{
			assert(partition_id_ != INVALID_INDEX);
			assert(n_partitions_ > 0);
			assert(global_id_.size() == owner_id_.size());
			assert(global_id_.size() == global_to_local_.size());

			if(partition_id_ == INVALID_INDEX) return false;
			if(n_partitions_ <= 0) return false;

			for(Integer i = 0; i < global_id_.size(); ++i) {
				assert(global_id_[i] != INVALID_INDEX);
				
				if(global_id_[i] == INVALID_INDEX) {
					return false;
				}
			}

			for(Integer i = 0; i < owner_id_.size(); ++i) {
				assert(owner_id_[i] != INVALID_INDEX);
				
				if(owner_id_[i] == INVALID_INDEX) {
					return false;
				}
			}

			std::vector<bool> visited(global_id_.size(), false);
			for(auto gtl : global_to_local_) {
				visited[gtl.second] = true;
				assert(gtl.first != INVALID_INDEX);

				if(gtl.first == INVALID_INDEX) {
					return false;
				}
			}

			for(auto v : visited) {
				assert(v);
				if(!v) return false;
			}

			if(!check_partitions) return true;

			assert(partitions_.size() == global_id_.size());

			for(Integer i = 0; i < partitions_.size(); ++i) {
				assert(partitions_[i].size() > 0);

				auto it = std::find(partitions_[i].begin(), partitions_[i].end(), partition_id_);
				assert(it != partitions_[i].end());

				if(it == partitions_[i].end()) {
					return false;
				}
			}

			return true;
		}

	private:
		Integer partition_id_;
		Integer n_partitions_;
		std::vector<Integer> global_id_;
		std::vector<Integer> owner_id_;

		//global to local
		std::map<Integer, Integer> global_to_local_;
		std::vector< std::vector<Integer> > partitions_;

		// std::vector<Integer> n_x_partition_;
	};

}

#endif //MARS_DOF_MAP_HPP
