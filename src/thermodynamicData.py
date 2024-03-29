import numpy as np
import scipy as sp

#=====================================================================
class ThermodynamicData:
#=====================================================================
    """ 
    La classe pour gérer les données thermochimiques des tables JANAF
    """
    
    deltah0 = 275.35e3
    deltas0 = 69.50


    R = 8.314


    Ttab = [298.15, 300, 350, 400, 450, 500, 
             600,  700,  800,  900, 1000,
            1100, 1200, 1300, 1400, 1500, 
            1600, 1700, 1800, 1900, 2000, 
            2100, 2200, 2300, 2400, 2500,
            2600, 2700, 2800, 2900, 3000,
            3100, 3200, 3300, 3400, 3500,
            3600, 3700, 3800, 3900, 4000,
            4100, 4200, 4300, 4400, 4500,
            4600, 4700, 4800, 4900, 5000,
            5100, 5200, 5300, 5400, 5500,
            5600, 5700, 5800, 5900, 6000]


        #--------------------------
        #           O
        #--------------------------
    
    
    ttab_o = [298.15,
       300,  350,  400,  450,  500, 
        600,  700,  800,  900, 1000,
        1100, 1200, 1300, 1400, 1500, 
        1600, 1700, 1800, 1900, 2000, 
        2100, 2200, 2300, 2400, 2500,
        2600, 2700, 2800, 2900, 3000,
        3100, 3200, 3300, 3400, 3500,
        3600, 3700, 3800, 3900, 4000,
        4100, 4200, 4300, 4400, 4500,
        4600, 4700, 4800, 4900, 5000,
        5100, 5200, 5300, 5400, 5500,
        5600, 5700, 5800, 5900, 6000 ]

    s0tab_o = [161.058,
        161.194, 164.551, 167.430, 169.953, 172.197,
        176.060, 179.310, 182.116, 184.585, 186.790,
        188.782, 190.599, 192.270, 193.816, 195.254,
        196.599, 197.862, 199.053, 200.179, 201.247,
        202.263, 203.232, 204.158, 205.045, 205.896,
        206.714, 207.502, 208.261, 208.995, 209.704,
        210.391, 211.057, 211.704, 212.332, 212.943,
        213.537, 214.117, 214.682, 215.234, 215.772,
        216.299, 216.814, 217.318, 217.812, 218.295,
        218.769, 219.234, 219.690, 220.138, 220.578,
        221.010, 221.435, 221.853, 222.264, 222.668,
        223.065, 223.457, 223.842, 224.222, 224.596 ]


    deltahftab_o =  [ 249.173,
        249.187, 249.537, 249.868, 250.180, 250.474,
        251.013, 251.494, 251.926, 252.320, 252.682,
        253.018, 253.332, 253.627, 253.906, 254.171,
        254.421, 254.659, 254.884, 255.097, 255.299,
        255.488, 255.667, 255.835, 255.992, 256.139,
        256.277, 256.405, 256.525, 256.637, 256.741,
        256.838, 256.929, 257.014, 257.094, 257.169,
        257.241, 257.309, 257.373, 257.436, 257.496,
        257.554, 257.611, 257.666, 257.720, 257.773,
        257.825, 257.876, 257.926, 257.974, 258.021,
        258.066, 258.110, 258.150, 258.189, 258.224,
        258.255, 258.282, 258.304, 258.321, 258.332 ]

        #--------------------------
        #           CO2
        #--------------------------
    
    
    ttab_co2 = [298.15,
        300,    400,     500, 
        600,  700,  800,  900, 1000,
        1100, 1200, 1300, 1400, 1500, 
        1600, 1700, 1800, 1900, 2000, 
        2100, 2200, 2300, 2400, 2500,
        2600, 2700, 2800, 2900, 3000,
        3100, 3200, 3300, 3400, 3500,
        3600, 3700, 3800, 3900, 4000,
        4100, 4200, 4300, 4400, 4500,
        4600, 4700, 4800, 4900, 5000,
        5100, 5200, 5300, 5400, 5500,
        5600, 5700, 5800, 5900, 6000 ]

    s0tab_co2 = [ 213.795,
        214.025, 225.314, 234.901,
        243.283, 250.750, 257.494, 263.645, 269.299,
        274.528, 279.390, 283.932, 288.191, 292.199,
        295.983, 299.566, 302.968, 306.205, 309.293,
        312.244, 315.070, 317.781, 320.385, 322.890,
        325.305, 327.634, 329.885, 332.061, 334.169,
        336.211, 338.192, 340.116, 341.986, 343.804,
        345.574, 347.299, 348.979, 350.619, 352.219,
        353.782, 355.310, 356.803, 358.264, 359.694,
        361.094, 362.466, 363.810, 365.128, 366.422,
        367.691, 368.937, 370.161, 371.364, 372.547,
        373.709, 374.853, 375.979, 377.087, 378.178 ]


    deltahftab_co2 = [ -393.522,
        -393.523, -393.583, -393.666,
        -393.803, -393.983, -394.188, -394.405, -394.623,
        -394.838, -395.050, -395.257, -395.462, -395.668,
        -395.876, -396.060, -396.311, -396.542, -396.784,
        -397.039, -397.309, -397.596, -397.900, -398.222,
        -398.562, -398.921, -399.299, -399.695, -400.111,
        -400.545, -400.998, -401.470, -401.960, -402.467,
        -402.991, -403.532, -404.089, -404.662, -405.251,
        -405.856, -406.475, -407.110, -407.760, -408.426,
        -409.106, -409.802, -410.514, -411.242, -411.986,
        -412.746, -413.522, -414.314, -415.123, -415.949,
        -416.794, -417.658, -418.541, -419.445, -420.372 ]

        #--------------------------
        #           CO
        #--------------------------

    
    ttab_co = [ 298.15,
        300,    400,     500, 
        600,  700,  800,  900, 1000,
        1100, 1200, 1300, 1400, 1500, 
        1600, 1700, 1800, 1900, 2000, 
        2100, 2200, 2300, 2400, 2500,
        2600, 2700, 2800, 2900, 3000,
        3100, 3200, 3300, 3400, 3500,
        3600, 3700, 3800, 3900, 4000,
        4100, 4200, 4300, 4400, 4500,
        4600, 4700, 4800, 4900, 5000,
        5100, 5200, 5300, 5400, 5500,
        5600, 5700, 5800, 5900, 6000 ]

    s0tab_co = [ 197.653,
        197.833, 206.238, 212.831,
        218.319, 223.066, 227.277, 231.074, 234.538,
        237.726, 240.679, 243.431, 246.006, 248.426,
        250.707, 252.865, 254.912, 256.859, 258.714,
        260.486, 262.182, 263.808, 265.369, 266.871,
        268.318, 269.713, 271.060, 272.362, 273.623,
        274.844, 276.029, 277.178, 278.295, 279.382,
        280.438, 281.468, 282.471, 283.449, 284.403,
        285.335, 286.245, 287.135, 288.005, 288.856,
        289.690, 290.506, 291.306, 292.090, 292.859,
        293.613, 294.354, 295.080, 295.794, 296.495,
        297.184, 297.862, 298.528, 299.184, 299.829 ]


    deltahftab_co = [-110.527,
        -110.516, -110.102, -110.103,
        -110.150, -110.469, -110.905, -111.418, -111.983,
        -112.586, -113.217, -113.870, -114.541, -115.229,
        -115.933, -116.651, -117.384, -118.133, -118.896,
        -119.675, -120.470, -121.282, -122.109, -122.953,
        -123.813, -124.689, -125.582, -126.490, -127.415,
        -128.356, -129.312, -130.283, -131.270, -132.271,
        -133.288, -134.318, -135.363, -136.422, -137.495,
        -138.581, -139.681, -140.794, -141.921, -143.062,
        -144.215, -145.383, -146.563, -147.758, -148.967,
        -150.190, -151.427, -152.679, -153.946, -155.228,
        -156.527, -157.841, -159.172, -160.521, -161.887 ]
    
    
        #--------------------------
        #           CH4
        #--------------------------
    
    
    ttab_ch4 = [ 298.15,
        300,  350,   400,  450,   500, 
        600,  700,  800,  900, 1000,
        1100, 1200, 1300, 1400, 1500, 
        1600, 1700, 1800, 1900, 2000, 
        2100, 2200, 2300, 2400, 2500,
        2600, 2700, 2800, 2900, 3000,
        3100, 3200, 3300, 3400, 3500,
        3600, 3700, 3800, 3900, 4000,
        4100, 4200, 4300, 4400, 4500,
        4600, 4700, 4800, 4900, 5000,
        5100, 5200, 5300, 5400, 5500,
        5600, 5700, 5800, 5900, 6000 ]

    s0tab_ch4 = [ 186.251,
        186.472, 192.131, 197.356, 202.291, 207.014,
        215.987, 224.461, 232.518, 240.205, 247.549,
        254.570, 261.287, 267.714,273.868, 279.763,
        285.413, 290.834, 296.039, 301.041, 305.853,
        310.485, 314.949, 319.255, 323.413, 327.431,
        331.317, 335.080, 338.725, 342.260, 345.690,
        349.021, 352.258, 355.406, 358.470, 361.453,
        364.360, 367.194, 369.959, 372.658, 375.293,
        377.868, 380.385, 382.846, 385.255, 387.612,
        389.921, 392.182, 394.399, 396.572, 398.703,
        400.794, 402.847, 404.861, 406.840, 408.784,
        410.694, 412.572, 414.418, 416.234, 418.020 ]

    deltahftab_ch4= [ -74.873,
        -74.929, -76.461, -77.969, -79.422, -80.802,
        -83.308, -85.452, -87.238, -88.692, -89.849,
        -90.750, -91.437, -91.945, -92.308, -92.553,
        -92.703, -92.780, -92.797, -92.770, -92.709,
        -92.624, -92.521, -92.409, -92.291, -92.174,
        -92.060, -91.954, -91.857, -91.773, -91.705,
        -91.653, -91.621, -91.609, -91.619, -91.654,
        -91.713, -91.798, -91.911, -92.051, -92.222,
        -92.422, -92.652, -92.914, -93.208, -93.533,
        -93.891, -94.281, -94.702, -95.156, -95.641,
        -96.157, -96.703, -97.278, -97.882, -98.513,
        -99.170, -99.852, -100.557, -101.284, -102.032 ]
    
        #--------------------------
        #           O2
        #--------------------------

    ttab_o2 = [ 298.15,
         300,  350,  400,  450,  500, 
         600,  700,  800,  900, 1000,
        1100, 1200, 1300, 1400, 1500, 
        1600, 1700, 1800, 1900, 2000, 
        2100, 2200, 2300, 2400, 2500,
        2600, 2700, 2800, 2900, 3000,
        3100, 3200, 3300, 3400, 3500,
        3600, 3700, 3800, 3900, 4000,
        4100, 4200, 4300, 4400, 4500,
        4600, 4700, 4800, 4900, 5000,
        5100, 5200, 5300, 5400, 5500,
        5600, 5700, 5800, 5900, 6000 ]

    s0tab_o2 = [ 205.147,
        205.329, 209.880, 213.871, 217.445, 220.693,
        226.451, 231.466, 235.921, 239.931, 243.378,
        246.922, 250.010, 252.878, 255.556, 258.068,
        260.434, 262.672, 264.796, 266.818, 268.748,
        270.595, 272.366, 274.069, 275.709, 277.290,
        278.819, 280.297, 281.729, 283.118, 284.466,
        285.776, 287.050, 288.291, 289.499, 290.677,
        291.826, 292.948, 294.044, 295.115, 296.162,
        297.186, 298.189, 299.171, 300.133, 301.076,
        302.002, 302.910, 303.801, 304.677, 305.538,
        306.385, 307.217, 308.037, 308.844, 309.639,
        310.424, 311.197, 311.960, 312.713, 313.457 ]

    deltahftab_o2 = [ 0,
        0, 0, 0, 0, 0,
        0, 0, 0, 0, 0,
        0, 0, 0, 0, 0,
        0, 0, 0, 0, 0,
        0, 0, 0, 0, 0,
        0, 0, 0, 0, 0,
        0, 0, 0, 0, 0,
        0, 0, 0, 0, 0,
        0, 0, 0, 0, 0,
        0, 0, 0, 0, 0,
        0, 0, 0, 0, 0,
        0, 0, 0, 0, 0 ]

        #--------------------------
        #           NH3
        #--------------------------

    ttab_nh3  = [ 298.15,
        300,        400,        500, 
        600,  700,  800,  900, 1000,
        1100, 1200, 1300, 1400, 1500, 
        1600, 1700, 1800, 1900, 2000, 
        2100, 2200, 2300, 2400, 2500,
        2600, 2700, 2800, 2900, 3000,
        3100, 3200, 3300, 3400, 3500,
        3600, 3700, 3800, 3900, 4000,
        4100, 4200, 4300, 4400, 4500,
        4600, 4700, 4800, 4900, 5000,
        5100, 5200, 5300, 5400, 5500,
        5600, 5700, 5800, 5900, 6000 ]

    s0tab_nh3 = [ 192.774,
        192.995, 203.663, 212.659,
        220.615, 227.829, 234.476, 240.669, 246.486,
        251.983, 257.199, 262.166, 266.907, 271.442,
        275.788, 279.957, 283.962, 287.815, 291.525,
        295.101, 298.552, 301.884, 305.104, 308.220,
        311.236, 314.158, 316.991, 319.740, 322.409,
        325.001, 327.521, 329.972, 332.358, 334.680,
        336.942, 339.147, 341.297, 343.395, 345.441,
        347.439, 349.391, 351.297, 353.161, 354.983,
        356.765, 358.508, 360.213, 361.882, 363.517,
        365.117, 366.685, 368.223, 369.732, 371.214,
        372.669, 374.098, 375.503, 376.883, 378.240 ]         



    deltahftab_nh3 = [ -45.898,
        -45.939, -48.041, -49.857,
        -51.374, -52.618, -53.621, -54.411, -55.013,
        -55.451, -55.746, -55.917, -55.982, -55.954,
        -55.847, -55.672, -55.439, -55.157, -54.833,
        -54.473, -54.084, -53.671, -53.238, -52.789,
        -52.329, -51.860, -51.386, -50.909, -50.433,
        -49.959, -49.491, -49.030, -48.578, -48.139,
        -47.713, -47.302, -46.908, -46.534, -46.180,
        -45.847, -45.539, -45.254, -44.996, -44.764,
        -44.561, -44.387, -44.242, -44.129, -44.047,
        -43.999, -43.979, -43.982, -44.006, -44.049,
        -44.112, -44.193, -44.291, -44.404, -44.531 ]


        #--------------------------
        #           N2
        #--------------------------

    ttab_n2 = [ 298.15,
        300,  350,  400,  450,  500, 
        600,  700,  800,  900, 1000,
        1100, 1200, 1300, 1400, 1500, 
        1600, 1700, 1800, 1900, 2000, 
        2100, 2200, 2300, 2400, 2500,
        2600, 2700, 2800, 2900, 3000,
        3100, 3200, 3300, 3400, 3500,
        3600, 3700, 3800, 3900, 4000,
        4100, 4200, 4300, 4400, 4500,
        4600, 4700, 4800, 4900, 5000,
        5100, 5200, 5300, 5400, 5500,
        5600, 5700, 5800, 5900, 6000 ]

    s0tab_n2 = [ 191.609,
        191.789, 196.281, 200.181, 203.633, 206.739,
        212.176, 216.866, 221.017, 224.757, 228.170,
        231.313, 234.226, 236.943, 239.487, 241.880,
        244.138, 246.275, 248.304, 250.234, 252.074,
        253.833, 255.517, 257.132, 258.684, 260.176,
        261.614, 263.001, 264.341, 265.637, 266.891,
        268.106, 269.285, 270.429, 271.541, 272.622,
        273.675, 274.699, 275.698, 276.671, 277.622,
        278.549, 279.456, 280.341, 281.208, 282.056, 
        282.885, 283.698, 284.494, 285.275, 286.041,
        286.792, 287.529, 288.253, 288.964, 289.662,
        290.348, 291.023, 291.687, 292.341, 292.984 ]        



    deltahftab_n2 = [ 0,
        0, 0, 0, 0, 0,
        0, 0, 0, 0, 0,
        0, 0, 0, 0, 0,
        0, 0, 0, 0, 0,
        0, 0, 0, 0, 0,
        0, 0, 0, 0, 0,
        0, 0, 0, 0, 0,
        0, 0, 0, 0, 0,
        0, 0, 0, 0, 0,
        0, 0, 0, 0, 0,
        0, 0, 0, 0, 0,
        0, 0, 0, 0, 0 ]

    
        #--------------------------
        #           H2
        #--------------------------

    
    ttab_h2 = [ 298.15,
        300,  350,  400,  450,  500, 
        600,  700,  800,  900, 1000,
        1100, 1200, 1300, 1400, 1500, 
        1600, 1700, 1800, 1900, 2000, 
        2100, 2200, 2300, 2400, 2500,
        2600, 2700, 2800, 2900, 3000,
        3100, 3200, 3300, 3400, 3500,
        3600, 3700, 3800, 3900, 4000,
        4100, 4200, 4300, 4400, 4500,
        4600, 4700, 4800, 4900, 5000,
        5100, 5200, 5300, 5400, 5500,
        5600, 5700, 5800, 5900, 6000 ]

    s0tab_h2 = [ 130.680,
        130.858, 135.325, 139.216, 142.656, 145.737,
        151.077, 155.606, 159.548, 163.051, 166.216,
        169.112, 171.790, 174.288, 176.633, 178.846,
        180.944, 182.940, 184.846, 186.669, 188.418,
        190.099, 191.718, 193.278, 194.785, 196.243,
        197.654, 199.021, 200.349, 201.638, 202.891,
        204.111, 205.299, 206.457, 207.587, 208.690,
        209.767, 210.821, 211.851, 212.860, 213.848,
        214.816, 215.765, 216.696, 217.610, 218.508,
        219.389, 220.255, 221.106, 221.943, 222.767,
        223.577, 224.374, 225.158, 225.931, 226.691,
        227.440, 228.177, 228.903, 229.619, 230.323 ]        
        

    deltahftab_h2 = [ 0,
        0, 0, 0, 0, 0,
        0, 0, 0, 0, 0,
        0, 0, 0, 0, 0,
        0, 0, 0, 0, 0,
        0, 0, 0, 0, 0,
        0, 0, 0, 0, 0,
        0, 0, 0, 0, 0,
        0, 0, 0, 0, 0,
        0, 0, 0, 0, 0,
        0, 0, 0, 0, 0,
        0, 0, 0, 0, 0,
        0, 0, 0, 0, 0 ]


        #--------------------------
        #           H2O
        #--------------------------
    
    ttab_h2o = [ 298.15,
        300,        400,        500, 
        600,  700,  800,  900, 1000,
        1100, 1200, 1300, 1400, 1500, 
        1600, 1700, 1800, 1900, 2000, 
        2100, 2200, 2300, 2400, 2500,
        2600, 2700, 2800, 2900, 3000,
        3100, 3200, 3300, 3400, 3500,
        3600, 3700, 3800, 3900, 4000,
        4100, 4200, 4300, 4400, 4500,
        4600, 4700, 4800, 4900, 5000,
        5100, 5200, 5300, 5400, 5500,
        5600, 5700, 5800, 5900, 6000     ]
    
    s0tab_h2o = [ 188.834,
        189.042, 198.788, 206.534,
        213.052, 218.739, 223.825, 228.459, 232.738,
        236.731, 240.485, 244.035, 247.407, 250.620,
        253.690, 256.630, 259.451, 262.161, 264.769,
        267.282, 269.706, 272.048, 274.312, 276.503,
        278.625, 280.683, 282.680, 284.619, 286.504,
        288.337, 290.120, 291.858, 293.550, 295.201,
        296.812, 298.384, 299.919, 301.420, 302.887, 
        304.322, 305.726, 307.101, 308.448, 309.767,
        311.061, 312.329, 313.574, 314.795, 315.993,
        317.171, 318.327, 319.464, 320.582, 321.682,
        322.764, 323.828, 324.877, 325.909, 326.926    ]
    
    deltahftab_h2o = [ -241.86,
        -241.844, -242.846, -243.826,
        -244.758, -245.632, -246.443, -247.185, -247.857,
        -248.460, -248.997, -249.473, -249.894, -250.265,
        -250.592, -250.881, -251.138, -251.368, -251.575,
        -251.762, -251.934, -252.092, -252.239, -252.379,
        -252.513, -252.643, -252.771, -252.897, -253.024,
        -253.152, -253.282, -253.416, -253.553, -253.696,
        -253.844, -253.997, -254.158, -254.326, -254.501,
        -254.684, -254.876, -255.078, -255.288, -255.508,
        -255.738, -255.978, -256.229, -256.491, -256.763,
        -257.046, -257.338, -257.639, -257.950, -258.268,
        -258.595, -258.930, -259.272, -259.621, -259.977 ]
    
    
    #-------------------------------------------------------------------------

#--------------------------------------------------------
    def __init__(self):
#--------------------------------------------------------

     
     self.data = {
    'h2' : { 't': self.ttab_h2  , 's0' :  self.s0tab_h2  , 'deltahf' : self.deltahftab_h2  },
    'n2' : { 't': self.ttab_n2  , 's0' :  self.s0tab_n2  , 'deltahf' : self.deltahftab_n2  },
    'h2o': { 't': self.ttab_h2o , 's0' :  self.s0tab_h2o , 'deltahf' : self.deltahftab_h2o },
    'nh3': { 't': self.ttab_nh3 , 's0' :  self.s0tab_nh3 , 'deltahf' : self.deltahftab_nh3 },
    'o2' : { 't': self.ttab_o2  , 's0' :  self.s0tab_o2  , 'deltahf' : self.deltahftab_o2  },
    'ch4': { 't': self.ttab_ch4 , 's0' :  self.s0tab_ch4 , 'deltahf' : self.deltahftab_ch4 },
    'co' : { 't': self.ttab_co  , 's0' :  self.s0tab_co  , 'deltahf' : self.deltahftab_co  },
    'co2': { 't': self.ttab_co2 , 's0' :  self.s0tab_co2 , 'deltahf' : self.deltahftab_co2 },
    'o'  : { 't': self.ttab_o   , 's0' :  self.s0tab_o   , 'deltahf' : self.deltahftab_o   }
     }


    def interpolation_lineaire(self, x, y, x0):
        try:
             f = sp.interpolate.interp1d(x, y, kind='linear')
             fx0 = f(x0)
        except:
            pass     
        return fx0

    def interpolation_lineaire(self, x, y, x0):
        try:
             f = sp.interpolate.interp1d(x, y, kind='linear')
             fx0 = f(x0)
        except ValueError as e:
             print("Erreur d'interpolation linéaire :", e)
             fx0 = None  
        except Exception as e:
             print("Une erreur s'est produite lors de l'interpolation linéaire :", e)
             raise
        else:
             return fx0
    
    def calculate(self,T):

        try:
             self.s0_h2      = self.interpolation_lineaire (self.ttab_h2 , self.s0tab_h2      , T )
             self.s0_n2      = self.interpolation_lineaire (self.ttab_n2 , self.s0tab_n2      , T )
             self.s0_h2o     = self.interpolation_lineaire (self.ttab_h2o, self.s0tab_h2o     , T )
             self.s0_nh3     = self.interpolation_lineaire (self.ttab_nh3, self.s0tab_nh3     , T )
             self.s0_o2      = self.interpolation_lineaire (self.ttab_o2,  self.s0tab_o2      , T )
             self.s0_ch4     = self.interpolation_lineaire (self.ttab_ch4, self.s0tab_ch4     , T )
             self.s0_co      = self.interpolation_lineaire (self.ttab_co,  self.s0tab_co      , T )
             self.s0_co2     = self.interpolation_lineaire (self.ttab_co2, self.s0tab_co2     , T )
             self.s0_o       = self.interpolation_lineaire (self.ttab_o  , self.s0tab_o       , T )

             self.deltah_h2  = self.interpolation_lineaire (self.ttab_h2 , self.deltahftab_h2 , T )
             self.deltah_n2  = self.interpolation_lineaire (self.ttab_n2 , self.deltahftab_n2 , T )
             self.deltah_h2o = self.interpolation_lineaire (self.ttab_h2o, self.deltahftab_h2o, T )
             self.deltah_nh3 = self.interpolation_lineaire (self.ttab_nh3, self.deltahftab_nh3, T )
             self.deltah_o2  = self.interpolation_lineaire (self.ttab_o2 , self.deltahftab_o2,  T )
             self.deltah_ch4 = self.interpolation_lineaire (self.ttab_ch4, self.deltahftab_ch4, T )
             self.deltah_co  = self.interpolation_lineaire (self.ttab_co , self.deltahftab_co,  T )
             self.deltah_co2 = self.interpolation_lineaire (self.ttab_co2, self.deltahftab_co2, T )
             self.deltah_o   = self.interpolation_lineaire (self.ttab_o  , self.deltahftab_o  , T )
        
        
             self.deltah2 = 0.5 * self.deltah_o2  +     self.deltah_co  -       self.deltah_co2 
             self.deltah3 =       self.deltah_co  + 2 * self.deltah_h2o -       self.deltah_ch4 - 1.5 * self.deltah_o2
             self.deltah4 =       self.deltah_h2o -     self.deltah_h2  - 0.5 * self.deltah_o2 
             self.deltah5 =       self.deltah_n2  + 3 * self.deltah_h2o - 2   * self.deltah_nh3 - 1.5 * self.deltah_o2
             self.deltah6 =   2 * self.deltah_o   -     self.deltah_o2 
             self.deltah7 =   2 * self.deltah_nh3 -     self.deltah_n2  - 3   * self.deltah_h2 
             self.deltah8 =       self.deltah_co  + 3 * self.deltah_h2  -       self.deltah_ch4 -       self.deltah_h2o 
             self.deltah9 =       self.deltah_co2 +     self.deltah_h2  -       self.deltah_co  -       self.deltah_h2o 

             self.deltah2 *= 1e3
             self.deltah3 *= 1e3
             self.deltah4 *= 1e3
             self.deltah5 *= 1e3
             self.deltah6 *= 1e3
             self.deltah7 *= 1e3
             self.deltah8 *= 1e3
             self.deltah9 *= 1e3

             self.deltas2 = 0.5 * self.s0_o2  +     self.s0_co  -       self.s0_co2 
             self.deltas3 =       self.s0_co  + 2 * self.s0_h2o -       self.s0_ch4 - 1.5 * self.s0_o2
             self.deltas4 =       self.s0_h2o -     self.s0_h2  - 0.5 * self.s0_o2 
             self.deltas5 =       self.s0_n2  + 3 * self.s0_h2o - 2   * self.s0_nh3 - 1.5 * self.s0_o2
             self.deltas6 =   2 * self.s0_o   -     self.s0_o2 
             self.deltas7 =   2 * self.s0_nh3 -     self.s0_n2  - 3   * self.s0_h2 
             self.deltas8 =       self.s0_co  + 3 * self.s0_h2  -       self.s0_ch4 -       self.s0_h2o 
             self.deltas9 =       self.s0_co2 +     self.s0_h2  -       self.s0_co  -       self.s0_h2o 

             self.K0 = np.exp(-self.deltah0 / (self.R * T) + self.deltas0 / self.R  )
             self.K2 = np.exp(-self.deltah2 / (self.R * T) + self.deltas2 / self.R  )
             self.K3 = np.exp(-self.deltah3 / (self.R * T) + self.deltas3 / self.R  )
             self.K4 = np.exp(-self.deltah4 / (self.R * T) + self.deltas4 / self.R  )
             self.K5 = np.exp(-self.deltah5 / (self.R * T) + self.deltas5 / self.R  )
             self.K6 = np.exp(-self.deltah6 / (self.R * T) + self.deltas6 / self.R  )
             self.K7 = np.exp(-self.deltah7 / (self.R * T) + self.deltas7 / self.R  )
             self.K8 = np.exp(-self.deltah8 / (self.R * T) + self.deltas8 / self.R  )
             self.K9 = np.exp(-self.deltah9 / (self.R * T) + self.deltas9 / self.R  )
        
             self.deltaG0_0 = self.deltah0 - T * self.deltas0
             self.deltaG0_2 = self.deltah2 - T * self.deltas2
             self.deltaG0_3 = self.deltah3 - T * self.deltas3
             self.deltaG0_4 = self.deltah4 - T * self.deltas4
             self.deltaG0_5 = self.deltah5 - T * self.deltas5
             self.deltaG0_6 = self.deltah6 - T * self.deltas6
             self.deltaG0_7 = self.deltah7 - T * self.deltas7
             self.deltaG0_8 = self.deltah8 - T * self.deltas8
             self.deltaG0_9 = self.deltah9 - T * self.deltas9

             self.deltaG0_h2  = 1000*self.deltah_h2   - T * self.s0_h2
             self.deltaG0_h2o = 1000*self.deltah_h2o  - T * self.s0_h2o
             self.deltaG0_o   = 1000*self.deltah_o    - T * self.s0_o
             self.deltaG0_o2  = 1000*self.deltah_o2   - T * self.s0_o2
             self.deltaG0_co  = 1000*self.deltah_co   - T * self.s0_co
             self.deltaG0_co2 = 1000*self.deltah_co2  - T * self.s0_co2
             self.deltaG0_ch4 = 1000*self.deltah_ch4  - T * self.s0_ch4
             self.deltaG0_n2  = 1000*self.deltah_n2   - T * self.s0_n2
             self.deltaG0_nh3 = 1000*self.deltah_nh3  - T * self.s0_nh3
        except ValueError as e:
             print("Erreur de calcul :", e)
        except Exception as e:
             print("Une erreur s'est produite lors des calculs des valeurs thermodynamiques issues des tables JANAF :\n", e)
             raise
        else:
            pass    
      
        
    def print(self):

        for tt in self.Ttab:
            self.calculate(tt)
            print(f"T {tt:.2f} \t S0(H2)  {self.s0_h2:12.2f}  \t deltaH(H2)  {self.deltah_h2 : 12.2f}")

        for tt in self.Ttab:
            self.calculate(tt)
            print(f"T {tt:.2f} \t S0(N2)  {self.s0_n2:12.2f}  \t deltaH(N2)  {self.deltah_n2 : 12.2f}")

        for tt in self.Ttab:
            self.calculate(tt)
            print(f"T {tt:.2f} \t S0(H2O) {self.s0_h2o:12.3f} \t deltaH(H2O) {self.deltah_h2o: 12.3f}")

        for tt in self.Ttab:
            self.calculate(tt)
            print(f"T {tt:.2f} \t S0(NH3) {self.s0_nh3:12.3f} \t deltaH(NH3) {self.deltah_nh3: 12.3f}")

        for tt in self.Ttab:
            self.calculate(tt)
            print(f"T {tt:.2f} \t S0(O2)  {self.s0_o2: 12.3f} \t deltaH(O2)  {self.deltah_o2 : 12.3f}")
        
        for tt in self.Ttab:
            self.calculate(tt)
            print(f"T {tt:.2f} \t S0(CH4) {self.s0_ch4:12.3f} \t deltaH(CH4) {self.deltah_ch4: 12.3f}")

        for tt in self.Ttab:
            self.calculate(tt)
            print(f"T {tt:.2f} \t S0(CO)  {self.s0_co: 12.3f} \t deltaH(CO)  {self.deltah_co : 12.3f}")
            
        for tt in self.Ttab:
            self.calculate(tt)
            print(f"T {tt:.2f} \t S0(CO2) {self.s0_co2:12.3f} \t deltaH(CO2) {self.deltah_co2: 12.3f}")

        for tt in self.Ttab:
            self.calculate(tt)
            print(f"T {tt:.2f} \t S0(O)   {self.s0_o :12.2f}  \t deltaH(O)   {self.deltah_o  : 12.2f}")
            
    def print7(self):
        for tt in self.Ttab:
            self.calculate(tt)
            print(f"T {tt:.2f} \t deltaH(7)   {self.deltah7 :12.2f}  \t deltaS(7)   {self.deltas7  : 12.2f}\t K7   {self.K7  : 12.2g}")

    def print8(self):
        for tt in self.Ttab:
            self.calculate(tt)
            print(f"T {tt:.2f} \t deltaH(8)   {self.deltah8 :12.2f}  \t deltaS(8)   {self.deltas8  : 12.2f}\t K8   {self.K8  : 12.2g}")

    def print9(self):
        for tt in self.Ttab:
            self.calculate(tt)
            print(f"T {tt:.2f} \t deltaH(9)   {self.deltah9 :12.2f}  \t deltaS(9)   {self.deltas9  : 12.2f}\t K8   {self.K9  : 12.2g}")

    def buildEllinghmaDiagram(self):
        file = open ("data.txt", "w")
        file.write("IGOR\n")
        file.write("WAVES Tk deltaG0 deltaG2 deltaG3 deltaG4 deltaG5 deltaG6 deltaG7 deltaG8 deltaG9\n")
        file.write("BEGIN\n")
        for tt in self.Ttab:
            self.calculate(tt)
            file.write(f"{tt:10.2f}\t{self.deltaG0_0:10.2f}\t{self.deltaG0_2:10.2f}\t{self.deltaG0_3:10.2f}\t{self.deltaG0_4:10.2f}\t")
            file.write(f"{self.deltaG0_5:10.2f}\t{self.deltaG0_6:10.2f}\t{self.deltaG0_7:10.2f}\t{self.deltaG0_8:10.3f}\t{self.deltaG0_9:10.2f}\n")
        file.write("END\n")

    def printJANAFData (self, t):
        self.calculate(t)
        print (f"H2   deltaH0 {self.deltah_h2:10.3f} kJ/mol S0 {self.s0_h2:10.3f} J/K/mol deltaG0 {self.deltaG0_h2/1000:10.3f} kJ/mol")
        print (f"H2O  deltaH0 {self.deltah_h2o:10.3f} kJ/mol S0 {self.s0_h2o:10.3f} J/K/mol deltaG0 {self.deltaG0_h2o/1000:10.3f} kJ/mol")
        print (f"O    deltaH0 {self.deltah_o:10.3f} kJ/mol S0 {self.s0_o:10.3f} J/K/mol deltaG0 {self.deltaG0_o/1000:10.3f} kJ/mol")
        print (f"O2   deltaH0 {self.deltah_o2:10.3f} kJ/mol S0 {self.s0_o2:10.3f} J/K/mol deltaG0 {self.deltaG0_o2/1000:10.3f} kJ/mol")
        print (f"CO   deltaH0 {self.deltah_co:10.3f} kJ/mol S0 {self.s0_co:10.3f} J/K/mol deltaG0 {self.deltaG0_co/1000:10.3f} kJ/mol")
        print (f"CO2  deltaH0 {self.deltah_co2:10.3f} kJ/mol S0 {self.s0_co2:10.3f} J/K/mol deltaG0 {self.deltaG0_co2/1000:10.3f} kJ/mol")
        print (f"CH4  deltaH0 {self.deltah_ch4:10.3f} kJ/mol S0 {self.s0_ch4:10.3f} J/K/mol deltaG0 {self.deltaG0_ch4/1000:10.3f} kJ/mol")
        print (f"N2   deltaH0 {self.deltah_n2:10.3f} kJ/mol S0 {self.s0_n2:10.3f} J/K/mol deltaG0 {self.deltaG0_n2/1000:10.3f} kJ/mol")
        print (f"NH3  deltaH0 {self.deltah_nh3:10.3f} kJ/mol S0 {self.s0_nh3:10.3f} J/K/mol deltaG0 {self.deltaG0_nh3/1000:10.3f} kJ/mol")


    def printReaction (self, species):

        reaction = ""
        for key, value in species.items():
            if value < 0:
                if value == -1:
                    reaction += f"{key} + "
                else:
                    reaction += f"{-value} {key} + "

        reaction = reaction[:-3] + " -> "

        for key, value in species.items():
            if value > 0:
                if value == 1:
                    reaction += f"{key} + "
                else:
                    reaction += f"{value} {key} + "

        reaction = reaction[:-3]

        print()
        print(reaction)
        print()

    def calculateParametersForReaction (self,species,T):
        deltaG0 = 0
        deltaH0 = 0
        deltaS0 = 0

        ngas = 0

        for key, value in species.items():
            s0 = self.interpolation_lineaire (self.data[key]['t'] , self.data[key]['s0'], T )
            h0 = self.interpolation_lineaire (self.data[key]['t'] , self.data[key]['deltahf'], T )
            h0 *= 1000
            g0 = h0 - T * s0
            deltaH0 += h0 * value
            deltaS0 += s0 * value
            deltaG0 += h0 * value - T
            print (f"{key:5s}", f"h0 {h0/1000:10.3f} kJ/mol", f"   s0 {s0:8.3f} J/K/mol", f"   g0 {g0/1000:10.3f} kJ/mol")
            ngas = ngas + value

        T0 = deltaG0 / deltaS0
        logK = -deltaG0/8.31/T
        K = np.exp(logK)

        print()

        print(f"deltaH0 : {deltaH0/1000:10.3f} kJ")
        print(f"deltaS0 : {deltaS0:10.3f}")
        print(f"deltaG0 : {deltaG0/1000:10.3f} kJ")
        print(f"T0      : {T0:10.3f} K")
        print(f"logK    : {logK:10.3g}")
        print(f"K       : {K:10.3g} bar^{ngas:4.2f}")


        

        
