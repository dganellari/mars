#ifndef MARS_COMMUNICATOR_H
#define MARS_COMMUNICATOR_H

#include "base.hpp"

#ifdef WITH_MPI

#include <mpi.h>
#include <assert.h>
#include <vector>
#include <string>
#include <ostream>



namespace mars {

	using byte = unsigned char;
	
	std::string mpi_error_2_string(const int errorCode);
	int world_rank();
	
	MPI_Datatype mpi_data_type_double();
	MPI_Datatype mpi_data_type_long();
	MPI_Datatype mpi_data_type_float();
	MPI_Datatype mpi_data_type_int();
	MPI_Datatype mpi_data_type_char();
	MPI_Datatype mpi_data_type_byte();
	
	///General template for mpi datatypes
	template< typename T >
	static MPI_Datatype mpi_data_type();
	
	template<>
	inline MPI_Datatype mpi_data_type<double>()
	{
		return mpi_data_type_double();
	}
	
	template<>
	inline MPI_Datatype mpi_data_type<long>()
	{
		return mpi_data_type_long();
	}
	
	template<>
	inline MPI_Datatype mpi_data_type<float>()
	{
		return mpi_data_type_float();
	}
	
	template<>
	inline MPI_Datatype mpi_data_type<int>()
	{
		return mpi_data_type_int();
	}
	
	template<>
	inline MPI_Datatype mpi_data_type<char>()
	{
		return mpi_data_type_char();
	}
	
	template<>
	inline MPI_Datatype mpi_data_type<byte>()
	{
		return mpi_data_type_byte();
	}
	
	class MPIOp {
	public:
		virtual ~MPIOp() {}
		
		inline MPI_Op value() const { return value_; }
		
	protected:
		explicit MPIOp(const MPI_Op &value);
		MPI_Op value_;
	};
	
	class MPIMax : public MPIOp {
	public:
		MPIMax();
	};
	
	class MPIMin : public MPIOp {
	public:
		MPIMin();
	};
	
	class MPISum : public MPIOp {
	public:
		MPISum();
	};
	
	class MPIProd : public MPIOp {
	public:
		MPIProd();
	};
	
	/**
		* \brief MPI wrapper for the communication
		*/
	class Communicator //: public Describable, public VerboseObject 
	{
	private:
		MPI_Comm mpi_comm_;
		int rank_, num_procs_;
		std::vector<MPI_Request>  send_reqs_, recv_reqs_;
		std::vector<int> destination_ranks_;
		bool verbose_;
		
	public:
		
		Communicator(const MPI_Comm mpi_comm);
		Communicator();
		
		Communicator(const Communicator &other)
		: mpi_comm_(other.mpi_comm_),
		rank_(other.rank_),
		num_procs_(other.num_procs_),
		send_reqs_(other.send_reqs_),
		recv_reqs_(other.recv_reqs_),
		destination_ranks_(other.destination_ranks_), //FIXME it is a bit dangerous
		verbose_(other.verbose_) 
		{ }
		
		~Communicator() {}
		
		inline bool verbose() const
		{
			return verbose_;
		}

		inline void set_verbose(const bool val) 
		{
			verbose_ = val;
		}

		void clear()
		{
			send_reqs_.clear();
			destination_ranks_.clear();
			recv_reqs_.clear();
		}
		
		Communicator & operator = (const Communicator &other)
		{
			if (this == &other) {
				return *this;
			}
			
			assert(!other.has_pending_requests());
			
			mpi_comm_ = other.mpi_comm_;
			send_reqs_ = other.send_reqs_;
			destination_ranks_ = other.destination_ranks_;
			recv_reqs_ = other.recv_reqs_;
			rank_ = other.rank_;
			num_procs_ = other.num_procs_;
			
			return *this;
		}
		
		/**
		 * @param buffer the outgoing message
		 * @param size the size of the outgoing message
		 * @param to the rank of the remote process
		 * @param tag the identifier of the communication
		 */
		template< typename T >
		int i_send(T * buffer, const int size, const int to, const int tag)
		{
			return i_send(buffer, size, to, tag, mpi_data_type<T>());
		}
		
		int i_send(void * buffer, const int size, const int to, const int tag, MPI_Datatype type);
		
		// template<class OutputStream>
		// void i_send(OutputStream &os, const int to, const int tag)
		// {
		// 	assert(os.size() > 0);
		// 	i_send(os.pointer(), static_cast<int> ( os.size() ), to, tag);
		// }
		
		/**
		 * @param buffer[OUT] the incoming message
		 * @param size the size of the incoming message
		 * @param to the rank of the remote process
		 * @param tag the identifier of the communication
		 */
		template< typename T >
		int i_recv(T * buffer, const int size, const int from, const int tag)
		{
			return i_recv(buffer, size, from, tag, mpi_data_type<T>());
		}
		
		int i_recv(void *buffer, const int size, const int from, const int tag, MPI_Datatype type);
		
		void exscan(void *input, void *output, const int n, MPI_Datatype type, MPI_Op op);
		
		template<typename T>
		void v_cum_sum(T *input, T *output, const int size)
		{
			T local_cum_sum = 0;
			for(int i = 0; i < size; ++i) local_cum_sum += input[i];
			
			T temp = 0;
			exscan((void *) &local_cum_sum, (void *) &temp, 1, mpi_data_type<T>(), MPISum().value());
			
			if(!size) return;
			
			output[0] = temp + input[0];
			for(int i = 1; i < size; ++i) {
				output[i] = output[i-1] + input[i];
			}
		}
		
		template<typename T>
		void exscan(T *input, T *output, const int size, MPI_Op op)
		{
			exscan(input, output, size, mpi_data_type<T>(), op);
		}
		
		template<typename T>
		void exscan(T *input, T *output, const int size, const MPIOp &op)
		{
			exscan(input, output, size, mpi_data_type<T>(), op.value());
		}
		
		template<typename T>
		void scan(T *input, T *output, const int size, const MPIOp &op)
		{
			scan(input, output, size, mpi_data_type<T>(), op.value());
		}
		
		template< typename T >
		void recv(T * buffer, const int size, const int from, const int tag)
		{
			recv(buffer, size, from, tag, mpi_data_type<T>());
		}
		
		///wait for all the communication to be finished
		inline void wait_all()
		{
			wait_all_recv();
			wait_all_send();
		}
		
		///@return true if not all the requests are completed
		inline bool has_pending_requests() const
		{
			return !recv_reqs_.empty() || !send_reqs_.empty();
		}
		
		///Waits for all the incoming communication to be finished
		void wait_all_recv();
		
		//Waits for all the outgoing communication to be finished
		void wait_all_send();
		
		/**
		 * @param send_buffer the data sent to all the other processes of this communicator
		 * @param recv_buffer the data received from all the other processes of this communicator
		 */
		template< typename T >
		void all_gather(T * send_buffer, T * recv_buffer, const int size) const
		{
			all_gather(send_buffer, recv_buffer, size, mpi_data_type<T>());
		}
		
		template< typename T >
		void broadcast(T * buffer, const int size, const int root) const
		{
			broadcast(buffer, size, root, mpi_data_type<T>() );
		}
		
		void all_reduce_in_place(void *values, const int size, MPI_Datatype type, MPI_Op op) const;
		void recv(void *buffer, const int size, const int from, const int tag, MPI_Datatype type);
		void scan(void *input, void *output, const int size, MPI_Datatype type, MPI_Op op);
		void all_gather(void * send_buffer, void * recv_buffer, const int size, MPI_Datatype type) const;
		void broadcast(void * buffer, const int size, const int root, MPI_Datatype type) const;
		
		/**
		 * @param the param[IN,OUT] and result of the global reduction
		 * @param op the associative operation performed
		 */
		template< typename T >
		void all_reduce(T * val, const int size, MPI_Op op) const
		{
			assert(size >= 0);
			all_reduce_in_place(val, size, mpi_data_type<T>(), op);
		}
		
		template< typename T >
		void all_reduce(T * val, const int size, const MPIOp &op) const
		{
			assert(size >= 0);
			all_reduce_in_place(val, size, mpi_data_type<T>(), op.value());
		}
		
		///@return true if the communicator has size 1
		inline bool is_alone() const
		{
			return size() == 1;
		}
		
		///@return the rank of this process in this communicator
		inline int rank() const
		{
			return rank_;
		}
		
		///@return the size of the communicator
		inline int size() const
		{
			return num_procs_;
		}
		
		///@return true if the rank of this process in this communicator is 0
		inline  bool is_root() const
		{
			return rank_ == 0;
		}
		
		///@return the wrapped mpi communicator object
		inline MPI_Comm get_mpi_comm() const
		{
			return mpi_comm_;
		}
		
		void set_mpi_comm(MPI_Comm mpi_comm_);
		
		/**
		 *	@param source the rank of the process for which we test if a send has been posted
		 *	@return the count of the data associated with size T
		 */
		template< typename T >
		int i_probe(const Integer source)
		{
			return i_probe(source, mpi_data_type<T>());
		}
		
		int i_probe(const Integer source, MPI_Datatype type);
		
		///@brief returns the number of incoming recv (use either wait_all_recv or i_probe to ensure recvs are completed)
		template<class OutputBuffer, class InputBuffer>
		Integer i_send_recv_all(
								std::vector<OutputBuffer> &send_buffer,
								std::vector<InputBuffer>  &recv_buffer)
		{
			std::vector<Integer> connected(size(), 0);
			
			const Integer n_buffers = send_buffer.size();
			for(Integer r = 0; r < n_buffers; ++r) {
				assert(send_buffer[r].empty() || r != rank());
				if(r == rank()) continue;
				connected[r] = !send_buffer[r].empty();
			}
			
			const Integer n_incoming = incoming(connected);
			
			for(Integer i = 0; i < size(); ++i) {
				if(send_buffer[i].empty()) continue;
				
				i_send(send_buffer[i], i, i);
			}
			
			for(Integer i = 0; i < n_incoming; ++i) {
				Integer rank, size;
				while ( !i_probe_any<byte>( &rank, &size ) ) {}
				recv_buffer[rank].reserve(size);
				i_recv(recv_buffer[rank].raw_ptr(), size, rank, this->rank());
			}
			
			return n_incoming;
		}
		
		void barrier() const;
		
		template<typename T>
		bool is_globally_exactly_equal(const T value)
		{
			T minmax[2] = { -value, value };
			all_reduce(minmax, 2, MPIMax());
			return (minmax[0] + minmax[1]) == 0;
		}
		
		/**
		 *	@tparam T the datatype which is expected from the message
		 *	@param source[OUT] the rank of the process that posted a message send
		 *	@return the count of the data associated with size T
		 */
		template< typename T >
		int i_probe_any(Integer * source)
		{
			return i_probe_any(source, mpi_data_type<T>());	
		}
		
		int i_probe_any(Integer * source, MPI_Datatype type);
		
		template< typename T >
		bool i_probe_any_with_tag(Integer * source, Integer * count, const Integer tag)
		{
			return i_probe_any_with_tag(source, count, tag, mpi_data_type<T>());
		}
		
		bool i_probe_any_with_tag(Integer * source, Integer * count, const Integer tag, MPI_Datatype type);
		
		
		template< typename T >
		bool i_probe_any(Integer * source, Integer * count, Integer * tag = NULL)
		{
			return i_probe_any(source,  count, tag, mpi_data_type<T>());
		}
		
		bool i_probe_any(Integer * source, Integer * count, Integer * tag, MPI_Datatype type);
		
		bool test_send_any(Integer * destination, Integer * index);
		bool test_recv_any(Integer * source, Integer * index);
		
		
		struct Activity {
			static const int SEND_COMPLETE = 0;
			static const int RECV_COMPLETE = 1;
			static const int PROBE = 2;
			static const int NONE = 3;
			
			Integer rank, index, size, tag;
			int type;
			
			Activity();
			void clear();
			
			inline bool empty() const
			{
				return rank < 0;
			}
			
			void describe(std::ostream &os) const;

			friend std::ostream &operator<<(std::ostream &os, const Activity &a)
			{
				a.describe(os);
				return os;
			}
		};
		
		///!Contains while true (something must happen otherwise it will go in deadlock)
		template<typename T>
		Activity test_activity() {
			Activity a;
			if(test_recv_any(&a.rank, &a.index)) {
				a.type = Activity::RECV_COMPLETE;
				assert(!a.empty());
				return a;
			}
			
			if(test_send_any(&a.rank, &a.index)) {
				a.type = Activity::SEND_COMPLETE;
				assert(!a.empty());
				return a;
			}
			
			if(i_probe_any<T>(&a.rank, &a.size, &a.tag)) {
				a.type = Activity::PROBE;
				assert(!a.empty());
				return a;
			}
			
			a.type = Activity::NONE;
			return a;
		}
		
		inline bool is_valid(const Integer rank)
		{
			return (rank >= 0) && (rank < size());
		}
		
		static bool is_null(Integer rank);
		
		static inline bool probe_null(const Integer count)
		{
			return count == -1;
		}
		
		void describe(std::ostream &os) const
		{
			os << "[" << rank() << "] ";
		}

		friend std::ostream &operator<<(std::ostream &os, const Communicator &comm)
		{
			comm.describe(os);
			return os;
		}
		
		void abort(const int errorCode = 0) const;
		
		inline Integer incoming(std::vector<Integer> &outgoing)
		{
			assert(outgoing[rank_] == 0);
			assert(Integer(outgoing.size()) == size());
			all_reduce(&outgoing[0], outgoing.size(), MPISum().value());
			return outgoing[rank()];
		}

		inline static std::ostream &logger();
		
		// bool unstructured_all_gather(ByteOutputStream &send_buffer, std::vector<ByteInputBuffer> &recv_buffer, const bool blocking);
	};

	template<class Fun>
	void serial_apply(
		const Communicator &comm,
		Fun fun)
	{
		for(Integer i = 0; i < comm.size(); ++i) {
			comm.barrier();
			if(i == comm.rank()) {
				fun();
			}
		} 

		comm.barrier();
	}


	template<class Fun>
	void root_apply(
		const Communicator &comm,
		Fun fun)
	{
		comm.barrier();
		
		if(comm.is_root()) {
			fun();
		}
		
		comm.barrier();
	}
}

#endif //MARS_COMMUNICATOR_H
#endif //WITH_MPI