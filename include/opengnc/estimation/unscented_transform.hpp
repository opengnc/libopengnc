#ifndef OPENGNC_ESTIMATION_UNSCENTED_TRANSFORM_HPP
#define OPENGNC_ESTIMATION_UNSCENTED_TRANSFORM_HPP

#include <stdexcept>
#include <Eigen/Core>
#include <boost/asio.hpp>

namespace opengnc {
namespace estimation {

template <typename model_type>
class unscented_transform {
public:
    enum
    {
        nx = model_type::input_length,
        ny = model_type::output_length ,
        nsigma = 2*model_type::input_length +1
    };

    typedef typename model_type::input_scalar_type input_scalar_type;
    typedef typename model_type::input_scalar_type output_scalar_type;

    typedef Eigen::Matrix<input_scalar_type, nx, 1> x_vec;
    typedef Eigen::Matrix<input_scalar_type, nx, nx> x_mat;
    typedef Eigen::Matrix<output_scalar_type, ny, 1> y_vec;
    typedef Eigen::Matrix<output_scalar_type, ny, ny> y_mat;
    typedef Eigen::Matrix<output_scalar_type, nx, ny> xy_mat;
    typedef Eigen::Matrix<output_scalar_type, nsigma, 1> w_vec;
    typedef Eigen::Matrix<output_scalar_type, nsigma, nsigma> w_mat;
    typedef Eigen::Matrix<input_scalar_type, nx, nsigma> sigma_x_mat;
    typedef Eigen::Matrix<output_scalar_type, ny, nsigma> sigma_y_mat;

    unscented_transform()
     : _init(false)
     , _nsigma(nsigma)
     , _alpha(1.0)
     , _beta(2.0)
     , _gamma(0)
     , _kappa(0)
     , _lambda(0)
    { }

    unscented_transform(const x_vec& x, const x_mat& Px)
        : unscented_transform()
    {
        init(x,Px);
    }

    // Recommend 1e-4 to 1
    void set_alpha(double value) { _alpha = value; }

    void set_beta(double value) { _beta = value; }

    // 0 for state estimation, 3-n for parameter estimation
    void set_kappa(double value) { _kappa = value; }

    void init(const x_vec& x, const x_mat& Px)
    {
        using namespace Eigen;

        int n = nx;
        if (nx == Eigen::Dynamic) n = x.size();

        if (_last_n != n) {
            _last_n = n;
            _lambda = std::pow(_alpha, 2) * (n + _kappa) - n;
            _gamma = std::sqrt(n + _lambda);
            _nsigma = 2*n+1;
        }

        x_mat S = Px.llt().matrixL();

        for (int i=0; i<n; i++) {
            _sigma_x.col(i) = x + _gamma*S.col(i);
            _sigma_x.col(i+n) = x - _gamma*S.col(i);
        }

        _sigma_x.col(2*n) = x;

        // mean weights
        double wm = 1.0/(2*(n+_lambda));
        _wm = VectorXd::Constant(_nsigma, wm);
        _wm(2*n)=_lambda/(n+_lambda);

        // covariance weights
        _wc = _wm;
        _wc(2*n) += 1.0 - _alpha*_alpha + _beta;

        _init = true;
    }

    void propagate(model_type& model, y_vec& y, y_mat& Py, xy_mat& Pxy)
    {
        _propagate(model, y, Py, Pxy, true);
    }

    void propagate(model_type& model, y_vec& y, y_mat& Py)
    {
        xy_mat Pxy;
        _propagate(model, y, Py, Pxy, false);
    }

protected:
    void _propagate(model_type& model, y_vec& y, y_mat& Py, xy_mat& Pxy, bool calc_cross)
    {
        if (!_init) {
            throw std::runtime_error("Unscented transform must be initialised before calling propagate.");
        }

        if (nx == Eigen::Dynamic) {
            _sigma_y.resize(ny, _sigma_x.cols());
        }

        for (int i = 0; i <_sigma_x.cols(); ++i) {
          const x_vec& xi = static_cast<const x_vec&>(_sigma_x.col(i));
          y_vec yi = model(xi);

          if (yi.size() > _sigma_y.rows()) {
              _sigma_y.conservativeResize(yi.size(), xi.size());

              for (int col = 0; col < _sigma_y.cols(); ++col)
              for (int row = _sigma_y.rows(); row < yi.size(); ++row) {
                  _sigma_y(row, col) = 0;
              }
          }

          _sigma_y.block(0, i, yi.size(), 1) = yi;
        }


        y = _sigma_y * _wm;
        sigma_y_mat Delta_sigma_y = _sigma_y.colwise() - y;
        Py = Delta_sigma_y * _wc.asDiagonal() * Delta_sigma_y.transpose();

        if (calc_cross) {
            x_vec& x = static_cast<const x_vec&>(_sigma_x.col(_sigma_x.cols() - 1));
            sigma_x_mat Delta_sigma_x = _sigma_x.colwise() - x;
            Pxy = Delta_sigma_x * _wc.asDiagonal() * Delta_sigma_y.transpose();
        }
    }

    bool _init;
    int _last_n;
    int _nsigma;
    double _alpha;
    double _beta;
    double _gamma;
    double _kappa;
    double _lambda;
    sigma_x_mat _sigma_x;
    sigma_y_mat _sigma_y;
    w_vec _wm;
    w_vec _wc;
};

}
}
#endif // OPENGNC_ESTIMATION_UNSCENTED_TRANSFORM_HPP
