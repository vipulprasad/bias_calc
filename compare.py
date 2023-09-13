import numpy as np
import matplotlib.pyplot as plt
import scienceplots
plt.style.use('science')

input_path = "/home/vipul/vipul/halo_clutering/bias_calc/box7500/z3.0/bias_files/"
output_path = "/home/vipul/vipul/halo_clutering/bias_calc/box7500/z3.0/plots/"

def plot(file):

        data = np.loadtxt(file)
        plt.errorbar(data1[:,0], data1[:,1], fmt = '.', yerr = data1[:,2], label = '1')
        #plt.title('mass bin: {}, {}'.format(bin1, bin2))
        #plt.figsize(6,6)
        #plt.xscale('log')
        #plt.yscale('log')
        plt.ylim((1 , 4))
        plt.xlabel(r'$r$')
        plt.ylabel(r'$b(r)$')
        plt.legend()
        plt.savefig('bias_plots/bias_{}'.format(file.split('\n')[0].split('M')[1][0:-3]+'png'))
        plt.show()


def compare(file1):
        print(file1)

        bin1 = file1.split('_')[1]+'_'+file1.split('_')[2]
        bin2 = file1.split('_')[3]+'_'+file1.split('_')[4][0:-4]

        data1 = np.loadtxt(input_path+file1)
        data2 = np.loadtxt(input_path+'bias_{}_{}.txt'.format(bin1, bin1))
        data3 = np.loadtxt(input_path+'bias_{}_{}.txt'.format(bin2, bin2))

        ymin = min(data1[:,1])-1
        ymax = max(data1[:,2])+1
        cross_bias = np.sqrt(data2[:,1]*data3[:,1])

        plt.figure(figsize=(8, 6))
        #plt.scatter(distance, bias, s = 8)
        plt.errorbar(data1[:,0], data1[:,1], fmt = '.', yerr = data1[:,2], elinewidth = 2, linewidth = 8, label = r'$b(m_{1}, m_{2})$')
        plt.errorbar(data2[:,0], cross_bias, fmt = '.', yerr = data2[:,2], elinewidth = 2, linewidth = 8, label = r'$\sqrt{b(m_{1}) b(m_{2})}$')
        #plt.errorbar(data3[:,0], data3[:,1], fmt = '.', yerr = data3[:,2], label = '3')
        plt.title('mass bins: ({} - {}), ({} - {})'.format(bin1.split('_')[0], bin1.split('_')[1], bin2.split('_')[0], bin2.split('_')[1]))
        #plt.figsize(6,6)
        plt.xscale('log')
        #plt.yscale('log')
        #plt.ylim((ymin , ymax))
        plt.xlabel(r'$r$', fontsize = 'larger')
        plt.ylabel(r'$b(r)$', fontsize = 'larger')
        plt.legend()
        plt.savefig(output_path+'bias_compare_{}_{}.png'.format(bin1, bin2))
        #plt.show()
        plt.close()

file_names = open('bias_files_7500.txt', 'r').readlines()

for file_name in file_names:
        print(file_name)
        compare(file_name.split('\n')[0])

#compare('bias_3.99E+11_1.34E+12_7.38E+10_8.01E+10.txt')
