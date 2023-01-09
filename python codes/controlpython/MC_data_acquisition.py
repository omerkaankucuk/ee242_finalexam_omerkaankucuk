import mpcontrolPC

N = 4000
signal1, signal2 = mpcontrolPC.get_data_from_MC(N, 'com4', 1)
mpcontrolPC.save_to_file(signal1, file_name='pi_stepin.dat')
mpcontrolPC.save_to_file(signal2, file_name='pi_stepout.dat')