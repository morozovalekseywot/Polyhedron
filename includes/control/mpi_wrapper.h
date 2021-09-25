#ifndef AMR_PARALLEL_CONTROLLER_H
#define AMR_PARALLEL_CONTROLLER_H

#include <allstd.h>
#ifdef ENABLE_MPI
#include <mpi.h>

#if   SIZE_MAX == UCHAR_MAX
    #define MPI_SIZE_T MPI_UNSIGNED_CHAR

#elif SIZE_MAX == USHRT_MAX
    #define MPI_SIZE_T MPI_UNSIGNED_SHORT

#elif SIZE_MAX == UINT_MAX
    #define MPI_SIZE_T MPI_UNSIGNED

#elif SIZE_MAX == ULONG_MAX
    #define MPI_SIZE_T MPI_UNSIGNED_LONG

#elif SIZE_MAX == ULLONG_MAX
    #define MPI_SIZE_T MPI_UNSIGNED_LONG_LONG

#else
    #error "what is happening here?"
#endif
#endif


class ProcessesMesh;

class mpi {

public:
    mpi() = delete;
    mpi(const mpi& ex) = delete;
    mpi(mpi&& ex) = delete;

    static void init(int argc, char **argv);

    static void finilize();

    /**
     * Структура для хранения информации о том, сколько ячеек будет
     * передано какому процессу.
     */
    class CellsTransaction {
    public:
        CellsTransaction(int target, int cells_count);
        CellsTransaction(const CellsTransaction& transaction);

        int target() const;
        int cells_count() const;

    private:
        int m_target;         /// Целевой процесс
        int m_cells_count;    /// Количество ячеек к отправке
    };

    using CellsTransactions = vector<CellsTransaction>;

    static double min(double value);
    static size_t min(size_t value);
    static vector<double> min(const vector<double>& values);

    static double max(double value);
    static size_t max(size_t value);

    static double sum(double value);
    static size_t sum(size_t value);
    static vector<double> sum(const vector<double>& values);

    /** Собирает значния value со всех процессов в вектор */
    static vector<double> all_gather(double value);
    static vector<size_t> all_gather(size_t value);
    static vector<double> all_gather(const vector<double>& values);
    static vector<size_t> all_gather(const array<size_t, 3>& values);

    static void barrier();

    static int rank();
    static uint urank();
    static const string& srank();
    static int size();
    static uint usize();

    static bool is_master();

    class master_stream {
    public:
        template <class T>
        master_stream& operator<<(const T& val) {
            if (mpi::is_master())
                std::cout << val;
            return *this;
        }

        master_stream& operator<<(std::ostream& (*f)(std::ostream&)) {
            if (mpi::is_master())
                f(std::cout);
            return *this;
        }

        void flush() {
            if (mpi::is_master()) {
                std::cout.flush();
            }
        }
    };

    static master_stream cout;

private:
    static int g_rank;
    static int g_size;
};

#endif //AMR_PARALLEL_CONTROLLER_H
