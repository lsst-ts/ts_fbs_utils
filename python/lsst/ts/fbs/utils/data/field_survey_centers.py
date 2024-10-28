# Name, RA (deg), dec (deg)

targets = (
    ("Rubin_SV_095_-25", 95.0, -25.0),  # High stellar densty, low extinction
    ("Rubin_SV_125_-15", 125.0, -15.0),  # High stellar densty, low extinction
    ("DESI_SV3_R1", 179.60, 0.000),  # DESI, GAMA, HSC DR2, KiDS-N
    # (360 - 145., -25.),
    ("Rubin_SV_225_-40", 225.0, -40.0),  # 225 High stellar densty, low extinction
    ("DEEP_A0", 216, -12.5),  # DEEP Solar Systen
    ("Rubin_SV_250_2", 250.0, 2.0),  # 250 High stellar densty, low extinction
    # (219.80, -0.600), # DESI, GAMA, HSC DR2, KiDS-N
    # (360 - 53., -25),
    # (270.891667, -30.033889), # Baade's Window
    ("Rubin_SV_300_-41", 300.0, -41.0),  # High stellar densty, low extinction
    ("Rubin_SV_280_-48", 280.0, -48.0),  # High stellar densty, low extinction
    # ("DEEP_B1", 354.4, -1.0), # DEEP Solar System, B1j
    ("HSC-SSP Deep2-3", 352.360542, -0.844150),  # HSC SSP Deep2-3
    ("DEEP_B0", 310, -19),  # DEEP Solar System
    ("ELAIS_S1", 9.45, -44.0),  # ELAIS-S1 LSST DDF
    ("XMM_LSS", 35.708333, -4.75),  # LSST DDF
    ("ECDFS", 53.125, -28.1),  # ECDFS
    ("COSMOS", 150.1, 2.1819444444444445),  # COSMOS
    ("EDFS_A", 58.9, -49.315),  # EDFS_a
    ("EDFS_B", 63.6, -47.6),  # EDFS_b
    ("Abell 2744", 3.58, -30.58),  # JWST 50 sq. arcminutes to 5-sigma depth of 29
    ("Abell 2764", 5.713, -49.249),  # Euclid ERO target
    ("NGC 6822", 296.236, -14.788),  # Euclid ERO target
)


def get_sv_fields() -> dict:
    targets_dict = dict(
        zip([t[0] for t in targets], [{"RA": t[1], "Dec": t[2]} for t in targets])
    )
    return targets_dict
