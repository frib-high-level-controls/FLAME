#ifndef SCSI_LINEAR_H
#define SCSI_LINEAR_H

#include <ostream>
#include <math.h>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "base.h"

template<typename State>
struct LinearElementBase : public ElementVoid
{
    typedef State state_t;

    LinearElementBase(const Config& c)
        :ElementVoid(c)
        ,transfer(2,2)
    {}
    ~LinearElementBase() {}

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
        this->transfer = boost::numeric::ublas::identity_matrix<double>(2);
        this->transfer(0,1) = c.get<double>("length");
    }

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
               cos = ::cos(angle),
               sin = ::sin(angle);

        this->transfer(0,0) = this->transfer(1,1) = cos;
        this->transfer(0,1) = P*sin;
        this->transfer(1,0) = -sin/P;
    }

    virtual const char* type_name() const {return "dipole";}
};

template<typename State>
struct LinearThinQuad : LinearElementBase<State>
{
    typedef LinearElementBase<State> base_t;
    LinearThinQuad(const Config& c)
        :base_t(c)
    {
        double L = c.get<double>("length"),
               K = c.get<double>("strength", 1.0),
               sK, sKL, cos, sin;

        if(K<0.0) {
            // defocus
            K = -K;
            sK = sqrt(K);
            sKL = sK*L;
            cos = ::cosh(sKL);
            sin = ::sinh(sKL);
        } else {
            // focusing
            sK = sqrt(K);
            sKL = sK*L;
            cos = ::cos(sKL);
            sin = ::sin(sKL);
        }

        this->transfer(0,0) = this->transfer(1,1) = cos;
        this->transfer(0,1) = sin/sK;
        this->transfer(1,0) = sK*sin;
    }

    virtual const char* type_name() const {return "quad";}
};

#endif // SCSI_LINEAR_H
