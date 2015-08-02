#ifndef SCSI_LINEAR_H
#define SCSI_LINEAR_H

#include <ostream>
#include <math.h>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "base.h"

/** @brief An Element based on a simple Transfer matrix
 *
 * x' = M * x
 */
template<typename State>
struct LinearElementBase : public ElementVoid
{
    typedef State state_t;

    LinearElementBase(const Config& c)
        :ElementVoid(c)
        ,transfer(boost::numeric::ublas::identity_matrix<double>(6))
    {}
    virtual ~LinearElementBase() {}

    virtual void advance(StateBase& s) const
    {
        State& ST = static_cast<State&>(s);
        advanceT(ST);
    }

    virtual void show(std::ostream& strm) const
    {
        ElementVoid::show(strm);
        strm<<"Transfer: "<<transfer<<"\n";
    }

    typedef boost::numeric::ublas::matrix<double> value_t;

    value_t transfer;
private:
    void advanceT(State& s) const
    {
        using boost::numeric::ublas::prod;
        s.state = prod(transfer, s.state);
    }
};

template<typename State>
struct LinearDrift : LinearElementBase<State>
{
    typedef LinearElementBase<State> base_t;
    LinearDrift(const Config& c)
        :base_t(c)
    {
        this->transfer(State::L_Z, State::P_Z) = c.get<double>("length");
    }
    virtual ~LinearDrift() {}

    virtual const char* type_name() const {return "drift";}
};

template<typename State>
struct LinearThinDipole : LinearElementBase<State>
{
    typedef LinearElementBase<State> base_t;
    LinearThinDipole(const Config& c)
        :base_t(c)
    {
        double angle = c.get<double>("angle"), // in rad.
               P = c.get<double>("radius", 1.0),
               off = c.get<bool>("vertical", false) ? State::L_Y : State::L_X ,
               cos = ::cos(angle),
               sin = ::sin(angle);

        this->transfer(off,off) = this->transfer(off+1,off+1) = cos;
        this->transfer(off,off+1) = P*sin;
        this->transfer(off+1,off) = -sin/P;
    }
    virtual ~LinearThinDipole() {}

    virtual const char* type_name() const {return "dipole";}
};

template<typename State>
struct LinearThinQuad : LinearElementBase<State>
{
    typedef LinearElementBase<State> base_t;
    LinearThinQuad(const Config& c)
        :base_t(c)
    {
        double L  = c.get<double>("length"),
               K  = c.get<double>("strength", 1.0),
               sK = sqrt(K),
               sKL=sK*L,
               cos = ::cos(sKL),
               sin = ::sin(sKL),
               cosh = ::cosh(sKL),
               sinh = ::sinh(sKL);
        unsigned Fdir, Ddir;

        if(K<0.0) {
            // defocus in X, focus in Y
            Fdir = State::L_Y;
            Ddir = State::L_X;
        } else {
            // focus in X, defocus in Y
            Fdir = State::L_X;
            Ddir = State::L_Y;
        }

        this->transfer(Fdir,Fdir) = this->transfer(Fdir+1,Fdir+1) = cos;
        this->transfer(Fdir,Fdir+1) = sin/sK;
        this->transfer(Fdir+1,Fdir) = sK*sin;

        this->transfer(Ddir,Ddir) = this->transfer(Ddir+1,Ddir+1) = cosh;
        this->transfer(Ddir,Ddir+1) = sinh/sK;
        this->transfer(Ddir+1,Ddir) = sK*sinh;
    }
    virtual ~LinearThinQuad() {}

    virtual const char* type_name() const {return "quad";}
};

#endif // SCSI_LINEAR_H
