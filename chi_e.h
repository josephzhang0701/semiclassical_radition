struct ret_em_prob {double chi_r;double ax; double ay;double az;};
ret_em_prob Wr(const double* g, double t);
ret_em_prob Wr(
        const double* g,
        double t,
        double Ex, double Ey, double Ez,
        double Bx, double By, double Bz
);
//double IQED_div_Icl (double x);
//double Ispin_div_Icl (double x);

