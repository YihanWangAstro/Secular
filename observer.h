



#ifndef SECULAR_OBSERVER_
#define SECULAR_OBSERVER_

#include<fstream>

namespace secular{
  struct Stream_observer
  {
      Stream_observer(std::ostream &out, double dt) : dt_{dt}, t_out_{0.0}, f_out_{out}, switch_{secular::is_on(dt)} { }

      template<typename State>
      void operator()(State const&x , double t)
      {
          if(switch_ && t >= t_out_) {
                f_out_ << t << ' ' << x << "\r\n";
                t_out_ += dt_;
          }
      }
    private:
      double const dt_;
      double t_out_;
      std::ostream& f_out_;
      const bool switch_;
  };


  struct SMA_Determinator
  {
      SMA_Determinator(double a_coef, double a_min) : a_min_{a_min}, a_coef_{a_coef}, detect_{secular::is_on(a_min)}  { }

      template<typename State>
      bool operator()(State const&x , double t)
      {
          if(detect_){
              double a = secular::calc_a(a_coef_, x.L1x(), x.L1y(), x.L1z(), x.e1x(), x.e1y(), x.e1z());
              return a<= a_min_;
          } else {
              return false;
          }
      }
    private:
      double const a_min_;
      double const a_coef_;
      bool const detect_;
  };
}


#endif
