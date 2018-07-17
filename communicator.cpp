
#include "base.hpp"

#ifdef WITH_MPI

#include "communicator.hpp"
#include <iostream>

#define MARS_MPI_CATCH_ERROR(expr)  { int ret = (expr); if(MPI_SUCCESS != ret)  std::cerr << "[" << world_rank() << "][Error] mpi error " << mpi_error_2_string(ret) << " code: " << ret << "(" << __FILE__ << ":" << __LINE__ << ")" << std::endl; }

namespace mars {

	std::ostream &Communicator::logger()
	{
		return std::cout;
	}

	MPI_Datatype mpi_data_type_double()
	{
		return MPI_DOUBLE;
	}
	
	MPI_Datatype mpi_data_type_long()
	{
		return MPI_LONG;
	}
	
	MPI_Datatype mpi_data_type_float()
	{
		return MPI_FLOAT;
	}

	MPI_Datatype mpi_data_type_int()
	{
		return MPI_INT;
	}
	
	MPI_Datatype mpi_data_type_char()
	{
		return MPI_CHAR;
	}

	MPI_Datatype mpi_data_type_byte()
	{
		return MPI_BYTE;
	}

	static int N_UNSTRUCTURED_ALL_GATHER = 0;
	static int COMMUNICATOR_N_BARRIERS = 0;
	const static int COMMUNICATOR_CHECK_BARRIERS = 0;

	MPIOp::MPIOp(const MPI_Op &value) : value_(value) {}
	MPIMax::MPIMax() : MPIOp(MPI_MAX) {}
	MPIMin::MPIMin() : MPIOp(MPI_MIN) {}
	MPISum::MPISum() : MPIOp(MPI_SUM) {}
	MPIProd::MPIProd() : MPIOp(MPI_PROD) {}

	int world_rank()
	{
		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		return rank;
	}

	std::string mpi_error_2_string(const int error_code)
	{
		std::string str;
		str.resize(MPI_MAX_ERROR_STRING);
		int len = MPI_MAX_ERROR_STRING;
		MPI_Error_string(error_code, &str[0], &len);
		str.resize(len);
		return str;
	}

	Communicator::Communicator(const MPI_Comm mpi_comm)
	: mpi_comm_(mpi_comm), rank_(0), num_procs_(0), send_reqs_(), recv_reqs_(), destination_ranks_()
	{
		MARS_MPI_CATCH_ERROR( MPI_Comm_rank(mpi_comm_, &rank_) );
		MARS_MPI_CATCH_ERROR( MPI_Comm_size(mpi_comm_, &num_procs_) );

		assert(rank() < size());
	}

	Communicator::Communicator()
	: mpi_comm_(MPI_COMM_WORLD), rank_(0), num_procs_(0), send_reqs_(), recv_reqs_(), destination_ranks_()
	{
		MARS_MPI_CATCH_ERROR( MPI_Comm_rank(mpi_comm_, &rank_) );
		MARS_MPI_CATCH_ERROR( MPI_Comm_size(mpi_comm_, &num_procs_) );

		assert(rank() < size());
	}

	void Communicator::all_reduce_in_place(void *values, const int size, MPI_Datatype type, MPI_Op op) const
	{
		assert(size >= 0);
		MARS_MPI_CATCH_ERROR(MPI_Allreduce( MPI_IN_PLACE, values, size, type, op, mpi_comm_) );
	}

	void Communicator::exscan(void *input, void *output, const int n, MPI_Datatype type, MPI_Op op)
	{
		MPI_Exscan(input, output, n, type, op, mpi_comm_);
	}

	void Communicator::barrier() const
	{
		if(COMMUNICATOR_CHECK_BARRIERS) {

			if(is_root()) logger() << "X----------------------------------------------------X" << std::endl;
			MPI_Barrier(mpi_comm_);
			for(Integer i = 0; i < size(); ++i) {
				if(i == rank()) {
					logger() << *this << " n_barriers: " << COMMUNICATOR_N_BARRIERS << std::endl;
				}
				MPI_Barrier(mpi_comm_);
			}

			if(is_root()) logger() << "X----------------------------------------------------X" << std::endl;
		}

		MPI_Barrier(mpi_comm_);

		if(COMMUNICATOR_CHECK_BARRIERS) {
			++COMMUNICATOR_N_BARRIERS;

			if(COMMUNICATOR_N_BARRIERS > 21) {
				std::cerr << *this << " has more calls" << std::endl;
			}
		}
	}

	void Communicator::wait_all_recv()
	{
		if(recv_reqs_.empty()) return;

		if (recv_reqs_.size() == 1) {
			int flag;
			MARS_MPI_CATCH_ERROR(MPI_Test(&recv_reqs_[0], &flag, MPI_STATUS_IGNORE));
			if(!flag) {
				MARS_MPI_CATCH_ERROR(MPI_Wait(&recv_reqs_[0], MPI_STATUS_IGNORE));
			}
		} else {
			MARS_MPI_CATCH_ERROR(MPI_Waitall(static_cast< int >( recv_reqs_.size() ), &recv_reqs_[0], MPI_STATUSES_IGNORE));
		}

		recv_reqs_.clear();
	}

	//Waits for all the outgoing communication to be finished
	void Communicator::wait_all_send()
	{
		if(send_reqs_.empty()) return;

		if (send_reqs_.size() == 1) {
			int flag;
			MARS_MPI_CATCH_ERROR(MPI_Test(&send_reqs_[0], &flag, MPI_STATUS_IGNORE));

			if(!flag) {
				MARS_MPI_CATCH_ERROR(MPI_Wait(&send_reqs_[0], MPI_STATUS_IGNORE));
			}
		} else {

			MARS_MPI_CATCH_ERROR(MPI_Waitall(static_cast< int >( send_reqs_.size() ), &send_reqs_[0], MPI_STATUSES_IGNORE));
		}

		send_reqs_.clear();
		destination_ranks_.clear();
	}

	bool Communicator::test_send_any(Integer * destination, Integer * index)
	{
		if(send_reqs_.empty()) return false;

		*destination = MPI_PROC_NULL;
		int flag(0);
		int m_index;
		MPI_Status status;
		MARS_MPI_CATCH_ERROR( MPI_Testany( static_cast< int >( send_reqs_.size() ), &send_reqs_[0], &m_index, &flag, &status) );
		*index = m_index;
		if(!flag)
			return false;


		if(m_index != MPI_UNDEFINED) {
			assert(m_index < int(destination_ranks_.size()));
			*destination = destination_ranks_[m_index];
			return true;
		}

		return false;
	}

	bool Communicator::test_recv_any(Integer * source, Integer * index)
	{
		if(recv_reqs_.empty()) return false;

		*source = MPI_PROC_NULL;
		int flag(0);
		int m_index;
		MPI_Status status;
		MARS_MPI_CATCH_ERROR( MPI_Testany( static_cast< int >( recv_reqs_.size() ), &recv_reqs_[0], &m_index, &flag, &status) );
		*index = m_index;
		if(!flag)
			return false;

		*source = status.MPI_SOURCE;

		return *source>= 0;
	}

	// bool Communicator::unstructured_all_gather(ByteOutputStream &send_buffer, std::vector<ByteInputBuffer> &recv_buffer, const bool blocking)
	// {

	// 	assert(!has_pending_requests());

	// 	int tag = size() * 10 + N_UNSTRUCTURED_ALL_GATHER;
	// 	N_UNSTRUCTURED_ALL_GATHER = (N_UNSTRUCTURED_ALL_GATHER > 1e6)? 0 : (N_UNSTRUCTURED_ALL_GATHER + 1);

	// 	for(Integer i = 0; i < size(); ++i) {
	// 		if(i == rank())
	// 			continue;
	// 			// logger() << (*this) << " sent packet with size: " << send_buffer.size()  << " bytes to " << i << std::endl;
	// 		i_send(send_buffer.pointer(), send_buffer.size(), i, tag);
	// 	}

	// 	recv_buffer.resize(size(), ByteInputBuffer());
	// 	const Integer n_incoming = size()-1;
		
	// 	for (Integer index = 0; index < n_incoming; ++index) {
	// 		Integer remote_rank, size;
	// 		while ( !this->i_probe_any<byte>(&remote_rank, &size) ) {}
	// 			assert(size);

	// 		ByteInputStream &is = recv_buffer[remote_rank];
	// 		is.reserve(size);

	// 		// logger() << (*this) << " recv packet with size: " << size  << " from " << remote_rank << std::endl;
	// 		i_recv(is.pointer(), size, remote_rank, tag);
	// 	}        

	// 	if(blocking)	
	// 		wait_all();
	// 	return true;
	// }

	void Communicator::recv(void *buffer, const int size, const int from, const int tag, MPI_Datatype type)
	{
		MARS_MPI_CATCH_ERROR( MPI_Recv(buffer, size, type, from, tag, mpi_comm_, MPI_STATUS_IGNORE) );
	}

	void Communicator::scan(void *input, void *output, const int size, MPI_Datatype type, MPI_Op op)
	{
		MPI_Scan(input, output, size, type, op, mpi_comm_);
	}

	void Communicator::all_gather(void * send_buffer, void * recv_buffer, const int size, MPI_Datatype type) const
	{
		MARS_MPI_CATCH_ERROR(MPI_Allgather(send_buffer, size, type, recv_buffer, size, type, mpi_comm_));		
	}

	void Communicator::broadcast(void * buffer, const int size, const int root, MPI_Datatype type) const
	{
		MARS_MPI_CATCH_ERROR(MPI_Bcast( buffer, size, type, root, mpi_comm_ ));
	}

	int Communicator::i_recv(void *buffer, const int size, const int from, const int tag, MPI_Datatype type)
	{
		if(this->verbose()) {
			logger() << *this << " i_receving: " << size << " bytes from " << from << " (tag=" << tag << ")" << std::endl;
			logger() << std::flush;
		}

		assert(size > 0);
		const int index = recv_reqs_.size();
		recv_reqs_.push_back(MPI_REQUEST_NULL);
		MARS_MPI_CATCH_ERROR(MPI_Irecv(buffer, size, type, from, tag, mpi_comm_, &recv_reqs_.back()));
		return index;
	}

	int Communicator::i_send(void * buffer, const int size, const int to, const int tag, MPI_Datatype type)
	{
		if(this->verbose()) {
			logger() << *this << " i_sending: " << size << " bytes to " << to << " (tag=" << tag << ")" << std::endl;
			logger() << std::flush;
		}

		assert(size > 0);
		const int index = send_reqs_.size();
		send_reqs_.push_back(MPI_REQUEST_NULL);
		destination_ranks_.push_back(to);
		MARS_MPI_CATCH_ERROR(MPI_Isend(buffer, size, type, to, tag, mpi_comm_, &send_reqs_.back()));
		return index;
	}

	void Communicator::set_mpi_comm(MPI_Comm mpi_comm)
	{
		this->mpi_comm_ = mpi_comm;
		MARS_MPI_CATCH_ERROR( MPI_Comm_rank(mpi_comm_, &rank_) );
		MARS_MPI_CATCH_ERROR( MPI_Comm_size(mpi_comm_, &num_procs_) );
	}

	int Communicator::i_probe(const Integer source, MPI_Datatype type)
	{
		int flag(0);
		MPI_Status status;
		MARS_MPI_CATCH_ERROR( MPI_Iprobe( source, rank_, mpi_comm_, &flag, &status) );

		if(!flag)
			return -1;

		int count;
		MPI_Get_count(&status, type, &count);
		return count;
	}

	int Communicator::i_probe_any(Integer * source, MPI_Datatype type)
	{
		*source = MPI_PROC_NULL;

		int flag(0);
		MPI_Status status;
		MARS_MPI_CATCH_ERROR( MPI_Iprobe( MPI_ANY_SOURCE, rank_, mpi_comm_, &flag, &status) );

		if(!flag)
			return -1;

		*source = status.MPI_SOURCE;

		int count;
		MPI_Get_count(&status, type, &count);

		return count;

	}

	bool Communicator::i_probe_any_with_tag(Integer * source, Integer * count, const Integer tag, MPI_Datatype type)
	{
		*source = MPI_PROC_NULL;

		int flag(0);
		MPI_Status status;
		MARS_MPI_CATCH_ERROR( MPI_Iprobe( MPI_ANY_SOURCE, tag, mpi_comm_, &flag, &status) );

		if(!flag)
			return false;

		*source = status.MPI_SOURCE;

		int mCount;
		MPI_Get_count(&status, type, &mCount);
		*count = mCount;
		return true;
	}

	bool Communicator::i_probe_any(Integer * source, Integer * count, Integer * tag, MPI_Datatype type)
	{
		*source = MPI_PROC_NULL;

		int flag(0);
		MPI_Status status;
		MARS_MPI_CATCH_ERROR( MPI_Iprobe( MPI_ANY_SOURCE, MPI_ANY_TAG, mpi_comm_, &flag, &status) );

		if(!flag)
			return false;

		*source = status.MPI_SOURCE;
		// logger() << _rank << "<=" << *source << std::endl;

		if(tag) {
			*tag = status.MPI_TAG;
		}

		int mCount;
		MPI_Get_count(&status, type, &mCount);
		*count = mCount;
		return *source >= 0;
	}

	bool Communicator::is_null(Integer rank)
	{
		return rank == MPI_PROC_NULL;
	}

	void Communicator::abort(const int error_code) const {
		MPI_Abort(mpi_comm_, error_code);
	}

	Communicator::Activity::Activity()
	: rank(MPI_PROC_NULL), 
	index(-1), 
	size(-1), 
	tag(-1), 
	type(NONE)
	{}

	void Communicator::Activity::clear()
	{
		rank = MPI_PROC_NULL; 
		index = -1; 
		size = -1;
		tag = -1;
		type = NONE;
	}

	void Communicator::Activity::describe(std::ostream &os) const
	{
		static const char * typeNames[4] = { "SEND_COMPLETE",  "RECV_COMPLETE", "PROBE",  "NONE" };

		os << "Rank: " << rank << "\n";
		os << "Index: " << index << "\n";
		os << "Size: " << size << "\n";
		os << "Tag: " << tag << "\n";
		os << "Type: " << typeNames[type] << "=" << type << "\n";
	}
}

//clean up macros
#undef MARS_MPI_CATCH_ERROR
#endif //WITH_MPI