// Copyright (c) 2018 ERGO-Code. See license.txt for license.

#ifndef IPX_LP_SOLVER_H_
#define IPX_LP_SOLVER_H_

#include <memory>
#include "basis.h"
#include "control.h"
#include "ipm.h"
#include "iterate.h"
#include "model.h"
#include "presolver.h"
#include "simplex_iterate.h"
#include "solution.h"
#include "user_model.h"

namespace ipx {

class LpSolver {
public:
    LpSolver() : presolver_{user_model_, model_} {}

    // Loads an LP model in the form given in the reference documentation.
    // @num_var: number of variables, must be > 0.
    // @obj: size num_var array of objective coefficients.
    // @lb: size num_var array of variable lower bounds, can have -INFINITY.
    // @lb: size num_var array of variable upper bounds, can have +INFINITY.
    // @num_constr: number of constraints, must be >= 0.
    // @Ap, @Ai, @Ax: constraint matrix in CSC format; indices can be unsorted.
    // @rhs: size num_constr array of right-hand side entries.
    // @constr_type: size num_constr array of entries '>', '<' and '='.
    // Returns:
    //  0
    //  IPX_ERROR_argument_null
    //  IPX_ERROR_invalid_dimension
    //  IPX_ERROR_invalid_matrix
    //  IPX_ERROR_invalid_vector
    Int LoadModel(Int num_var, const double* obj, const double* lb,
                  const double* ub, Int num_constr, const Int* Ap,
                  const Int* Ai, const double* Ax, const double* rhs,
                  const char* constr_type);

    // Loads a primal-dual point as starting point for the IPM.
    // @x: size num_var array
    // @xl: size num_var array, must satisfy xl[j] >= 0 for all j and
    //      xl[j] == INFINITY iff lb[j] == -INFINITY.
    // @xu: size num_var array, must satisfy xu[j] >= 0 for all j and
    //      xu[j] == INFINITY iff ub[j] == INFINITY.
    // @slack: size num_constr array, must satisfy
    //         slack[i] == 0 if constr_type[i] == '='
    //         slack[i] >= 0 if constr_type[i] == '<'
    //         slack[i] <= 0 if constr_type[i] == '>'
    // @y: size num_constr array, must satisfy
    //     y[i] >= 0 if constr_type[i] == '>'
    //     y[i] <= 0 if constr_type[i] == '<'
    // @zl: size num_var array, must satsify zl[j] >= 0 for all j and
    //      zl[j] == 0 if lb[j] == -INFINITY
    // @zu: size num_var array, must satsify zu[j] >= 0 for all j and
    //      zu[j] == 0 if ub[j] == INFINITY
    // When a starting point was loading successfully (return value 0), then
    // the next call to Solve() will start the IPM from that point, except that
    // primal and dual slacks with value 0 are made positive if necessary. The
    // IPM will skip the initial iterations and start directly with basis
    // preconditioning.
    // At the moment using a starting point is not possible when the model is
    // dualized during preprocessing. See parameters to turn dualization off.
    // Returns:
    // 0                            success
    // IPX_ERROR_argument_null      an argument was NULL
    // IPX_ERROR_invalid_vector     a sign condition was violated
    Int LoadIPMStartingPoint(const double* x, const double* xl,
                             const double* xu, const double* slack,
                             const double* y, const double* zl,
                             const double* zu);

    // Solves the model that is currently loaded in the object.
    // Returns GetInfo().status.
    Int Solve();

    // Returns the solver info from the last call to Solve(). See the reference
    // documentation for the meaning of Info values.
    Info GetInfo() const;

    // Returns the final IPM iterate from the last call to Solve(). An iterate
    // is available if GetInfo().status_ipm != IPX_STATUS_not_run. If no iterate
    // is available, the method does nothing. In the first version an
    // InteriorSolution object is returned in @solution. In the second version
    // each of the pointer arguments must either be NULL or an array of
    // appropriate dimension. If NULL, the quantity is not returned. Returns -1
    // if no IPM iterate was available and 0 otherwise.
    Int GetInteriorSolution(InteriorSolution& solution) const;
    Int GetInteriorSolution(double* x, double* xl, double* xu, double* slack,
                            double* y, double* zl, double* zu) const;

    // Returns the basic solution and basis from the last call to Solve(). A
    // basic solution and basis are available if GetInfo().status_crossover ==
    // IPX_STATUS_optimal || GetInfo().status_crossover == IPX_STATUS_imprecise.
    // Otherwise the method does nothing. In the first version a BasicSolution
    // object is returned in @solution. In the second version each of the
    // pointer arguments must either be NULL or an array of appropriate
    // dimension. If NULL, the quantity is not returned. Returns -1 if no basic
    // solution was available and 0 otherwise.
    Int GetBasicSolution(BasicSolution& solution) const;
    Int GetBasicSolution(double* x, double* slack, double* y, double* z,
                         Int* cbasis, Int* vbasis) const;

    // Returns/sets all paramters. Without calling SetParameters(), the solver
    // uses the default values of a Parameters object.
    Parameters GetParameters() const;
    void SetParameters(Parameters new_parameters);

    // Reads/writes parameters to file. Returns 0 on success, -1 on failure.
    Int ReadParameters(const char* filename);
    Int WriteParameters(const char* filename);

    // Discards the model and solution (if any) but keeps the parameters.
    void ClearModel();

    // Discards the starting point (if any).
    void ClearIPMStartingPoint();

    // -------------------------------------------------------------------------
    // The remaining methods are for debugging.
    // -------------------------------------------------------------------------

    // Returns the current IPM iterate without postsolve. The method does
    // nothing when no iterate is available (i.e. when IPM was not started).
    // @x, @xl, @xu, @zl, @zu: either NULL or size num_cols_solver arrays.
    // @y: either NULL or size num_rows_solver array.
    // Returns -1 if no IPM iterate was available and 0 otherwise.
    Int GetIterate(double* x, double* y, double* zl, double* zu, double* xl,
                   double* xu);

    // Returns the current basis postsolved.
    // - If crossover terminated successfully, this is the basis returned by
    //   GetBasicSolution().
    // - If crossover failed, this is the basis at which failure occured.
    // - If crossover was not called, this is the basis from the IPM
    //   preconditioner.
    // - If no basis is available, the method does nothing.
    // @cbasis: either NULL or size num_constr array.
    // @vbasis: either NULL or size num_var array.
    // Returns -1 if no basis was available and 0 otherwise.
    Int GetBasis(Int* cbasis, Int* vbasis);

    // Returns the constraint matrix from the solver (including slack columns)
    // and the diagonal from the (1,1) block of the KKT matrix corresponding to
    // the current IPM iterate. The method does nothing when no IPM iterate is
    // available (i.e. when IPM was not started).
    // @AIp: either NULL or size num_cols_solver + 1 array.
    // @AIi, @AIx: either NULL or size num_entries_solver arrays.
    // (If any of the three arguments is NULL, the matrix is not returned.)
    // @g: either NULL or size num_cols_solver array.
    // Returns -1 if no IPM iterate was available and 0 otherwise.
    Int GetKKTMatrix(Int* AIp, Int* AIi, double* AIx, double* g);

    // (Efficiently) computes the number of nonzeros per row and column of the
    // symbolic inverse of the basis matrix.
    // @rowcounts, @colcounts: either NULL or size num_rows_solver arrays.
    // Returns -1 if no basis was available and 0 otherwise.
    Int SymbolicInvert(Int* rowcounts, Int* colcounts);

private:
    void ClearSolution();
    void Presolve();
    void InteriorPointSolve();
    void RunIPM();
    void LoadStartingPoint(IPM& ipm);
    void ComputeStartingPoint(IPM& ipm);
    void RunInitialIPM(IPM& ipm);
    void BuildStartingBasis();
    void RunMainIPM(IPM& ipm);
    void RunCrossover();
    void PrintSummary();

    Control control_;
    Info info_;
    UserModel user_model_;
    Model model_;
    Presolver presolver_;
    std::unique_ptr<Iterate> iterate_;
    std::unique_ptr<Basis> basis_;
    std::unique_ptr<SimplexIterate> simplex_iterate_;

    // Postsolved IPM solution and crossover solution (if available).
    std::unique_ptr<InteriorSolution> interior_solution_;
    std::unique_ptr<BasicSolution> basic_solution_;

    // IPM starting point provided by user.
    std::unique_ptr<InteriorSolution> ipm_start_;
};

}  // namespace ipx

#endif  // IPX_LP_SOLVER_H_
