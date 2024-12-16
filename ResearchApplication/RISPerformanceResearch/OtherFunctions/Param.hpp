#pragma once

struct Param
{
    double BSden;
    double RISden;
    double RISnumber;
    double radius_min;
    double radius_max;
    double UEs_cell;
    double noise_level{1.0e-13};
    double beam_probability {0.01};

    double Antenna_gain {6.332573977646111e-05};

    Param(double BSden_, double RISden_, double RISnumber_, double radius_min_, double radius_max_, double UE_cell_): \
        BSden(BSden_), RISden(RISden_), RISnumber(RISnumber_), radius_min(radius_min_), radius_max(radius_max_), UEs_cell(UE_cell_){}
};
