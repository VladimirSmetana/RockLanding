#include "airforce.h"

    airforce::airforce (double Mah)
    {

        CX = -0.65;
        CY = 0;
//        if ((Mah > 0) && (Mah <= 0.25)) CY = 2.8;
//        if ((Mah > 0.25) && (Mah <= 1.1)) CY = 2.8+0.447*(Mah-0.25);
//        if ((Mah > 1.1) && (Mah <= 1.6)) CY = 3.18 - 0.660 * (Mah - 1.1);
//        if ((Mah > 1.6) && (Mah <= 3.6)) CY = 2.85 + 0.350 * (Mah - 1.6);
//        if (Mah > 3.6) CY = 12.78/Mah;

//        if (Mah <= 0.8) { CX = 0.29; };
//        if (Mah > 0.8 && Mah < 1.068) { CX = Mah-0.51; };
//        if (Mah >= 1.068) { CX = 0.597/Mah; };
//        if (CX > 1) { CX = 1; };
    }

    double airforce::getCX()
    {
        return CX;
    }

    double airforce::getCY()
    {
        return CY;
    }
