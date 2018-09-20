from __future__ import division

import math
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

def read_file(file_name):

    distribution_names = []
    distribution_types = []
    distribution_entries = []

    file_handle = open(file_name, 'r')
    lines = file_handle.read().splitlines()
    
    for i in range(0, len(lines)):
        parsed_line = lines[i].split(' ')

        tmp_name = parsed_line[0].split('_')

        distribution_names.append(' '.join(tmp_name))
        triple_parsed = []
        sholl_append = []
        
        del parsed_line[0]
        
        distribution_types.append(parsed_line[0])
        del parsed_line[0]
        if parsed_line[0] != '':
            if parsed_line[0][0] == '[':
    #            del parsed_line[0]
    #            del parsed_line[-1]
                for entry in parsed_line:
                    if entry[0] == '[':
                        double_parsed = []
                        double_parsed.append(entry[1:-1])
                    elif entry[-1] == ']':
                        double_parsed.append(entry[:-1])
                        triple_parsed.append(double_parsed)
                    else:
                        double_parsed.append(entry[:-1])
                for j in range(0, len(triple_parsed)):
                    sholl_append.append([float(x) for x in triple_parsed[j]])
                distribution_entries.append(sholl_append)
            else:
                distribution_entries.append(map(float,parsed_line))
        else:
            distribution_entries.append([0])
    
    return distribution_names, distribution_types, distribution_entries

def plot_2_distributions(distribution_1, distribution_2, distribution_name):
    fig = plt.figure()
    bin_num = 25
    both_distributions = distribution_1 + distribution_2
    [hist_peaks, bin_edges] =np.histogram(both_distributions, bins=bin_num) #get the bin edges
    [hist_peaks1, bin_edges] =np.histogram(distribution_1, bins=bin_num)
    [hist_peaks2, bin_edges] =np.histogram(distribution_2, bins=bin_num)
    
    print bin_edges
    bin_edges = np.ndarray.tolist(bin_edges)
    
    del bin_edges[-1]
#    plt.hist(distribution_1, alpha=0.7, normed=True, bins=bin_edges, label=legend_name_1)
#    plt.hist(distribution_2, alpha=0.7, normed=True, bins=bin_edges, label=legend_name_2)
    plt.plot(bin_edges, hist_peaks1, color='black', linestyle='solid', label=legend_name_1)
    plt.plot(bin_edges, hist_peaks2, color='black', linestyle='dashed',label=legend_name_2)
    plt.xlabel(distribution_name, fontsize=18)
    plt.ylabel('Frequency', fontsize=18)

    xt = plt.xticks()[0]
    xmin, xmax = min(xt), max(xt)
    lnspc = np.linspace(xmin, xmax, len(distribution_2))
#    be, ce = stats.expon.fit(distribution_1)
#    pdf_expon = stats.expon.pdf(lnspc, be, ce)
#    plt.plot(lnspc, pdf_expon, label = "Exp")
#    m, s = stats.norm.fit(distribution_1)
    ag, bg, cg = stats.gamma.fit(distribution_2)
    pdf_gamma = stats.gamma.pdf(lnspc, ag, bg, cg)
#    pdf_norm = stats.norm.pdf(lnspc, m, s)
#    plt.plot(lnspc, pdf_gamma, label="Fitted Distribution", color='red')
#    plt.plot(lnspc, pdf_norm, label="Norm")
    
#    print("Exponential Coef", be, ce)
#    print("Normal Coef", m, s)
#    print("Gamma Coef", ag, bg, cg)


    distribution_1.sort()
    distribution_2.sort()
    
    W_stat, pvalue = stats.mannwhitneyu(distribution_1, distribution_2, alternative='two-sided')
    
    percentile1 = [0]*len(distribution_1)
    percentile2 = [0]*len(distribution_2)
    
    delta = 1
    max1 = math.ceil(distribution_1[-1])
    max2 = math.ceil(distribution_2[-1])
    index1 = np.linspace(0, max1, int(max1/delta)).tolist()
    index2 = np.linspace(0, max2, int(max2/delta)).tolist()
    cumulative1 = [0]*int(max1/delta)
    cumulative2 = [0]*int(max2/delta)
    
    for i in range(0, int(max1/delta)):
        matches = [ii for ii, x in enumerate(distribution_1) if x > (i+1)*delta]
#        print matches
        if matches == []:
            cumulative1[i] = 1
        else:
            cumulative1[i] = (matches[0]+1)/(len(distribution_1)+1)

    for i in range(0, int(max2/delta)):
        matches = [ii for ii, x in enumerate(distribution_2) if x > (i+1)*delta]
        if matches == []:
            cumulative2[i] = 1
        else:
            cumulative2[i] = (matches[0]+1)/(len(distribution_2)+1)

    for i in range(0, len(distribution_1)):
        percentile1[i] = (i+1)/(len(distribution_1)+1)
        
    for j in range(0, len(distribution_2)):
        percentile2[j] = (j+1)/(len(distribution_2)+1)


#    fig = plt.figure()
#    plt.plot(index1, cumulative1, color='red')
#    plt.plot(index2, cumulative2, color='blue')
#    plt.show()
#    
#    fig = plt.figure()
#    plt.plot(distribution_1, percentile1, color='red')
#    plt.plot(distribution_2, percentile2, color='blue')
##    plt.yscale('log')
##    plt.xscale('log')
#    plt.show()
    
    KS_stat, p_val = stats.ks_2samp(distribution_1, distribution_2)
    
    n1 = len(distribution_1)
    n2 = len(distribution_2)
    z_u = abs(W_stat - (n1*n2/2))/math.sqrt(n1*n2*(n1+n2+1)/12)
    
    print '\n', distribution_name, '\nW stat:', W_stat, '\n', 'p value:', pvalue, '\nz value', z_u


    #plot KS-plot
    
#
    plt.legend()
    plt.show()
    
def plot_3_distributions(distribution_1, distribution_2, distribution_3, distribution_name):
    
    both_distributions = distribution_1 + distribution_2 + distribution_3
    bin_edges=np.histogram(both_distributions, bins=50)[1] #get the bin edges

    
    plt.hist(distribution_1, alpha=0.7, normed=True, bins=bin_edges, label=legend_name_1)
    plt.hist(distribution_2, alpha=0.7, normed=True, bins=bin_edges, label=legend_name_2)
    plt.hist(distribution_3, alpha=0.7, normed=True, bins=bin_edges, label=legend_name_3)
    plt.xlabel(distribution_name)
    plt.ylabel('Normalized Frequency')
    

    plt.legend()
    plt.show()    
    
def plot_4_distributions(distribution_1, distribution_2, distribution_3, distribution_4, distribution_name):
    
    both_distributions = distribution_1 + distribution_2 + distribution_3 + distribution_4
    bin_edges=np.histogram(both_distributions, bins=50)[1] #get the bin edges

    
    plt.hist(distribution_1, alpha=0.5, normed=True, bins=bin_edges, label=legend_name_1)
    plt.hist(distribution_2, alpha=0.5, normed=True, bins=bin_edges, label=legend_name_2)
    plt.hist(distribution_3, alpha=0.5, normed=True, bins=bin_edges, label=legend_name_3)
    plt.hist(distribution_4, alpha=0.5, normed=True, bins=bin_edges, label=legend_name_4)
    plt.xlabel(distribution_name)
    plt.ylabel('Normalized Frequency')
    

    plt.legend()
    plt.show()        
    
def plot_conditional_relationships(distribution_names, distributions, basic_parameters, conditional_parameters):
    for i in range(0, len(basic_parameters)):
        for j in range(0, len(conditional_parameters)):
            basic_index = basic_parameters[i]
            conditional_index = conditional_parameters[j]
            if len(distributions[conditional_index]) == len(distributions[basic_index]):
                
                print (distribution_names[basic_index] + ' depending on ' + distribution_names[conditional_index])
                plt.scatter(distributions[conditional_index], distributions[basic_index])
                plt.show()

def compare_conditional_relationships(distribution_names, distributions_1, distributions_2, basic_parameters, conditional_parameters):
    for i in range(0, len(basic_parameters)):
        for j in range(0, len(conditional_parameters)):
            basic_index = basic_parameters[i]
            conditional_index = conditional_parameters[j]
            if len(distributions_1[conditional_index]) == len(distributions_1[basic_index]):
                
                print (distribution_names[basic_index] + ' depending on ' + distribution_names[conditional_index])
                plt.scatter(distributions_1[conditional_index], distributions_1[basic_index], alpha=0.6)
                plt.scatter(distributions_2[conditional_index], distributions_2[basic_index], alpha=0.6)
                plt.show()

file_name_1 = '/home/zzchou/Documents/Codes/claiborne_distributions.txt'
#file_name_4 = '/home/zzchou/Documents/Codes/beining_distributions.txt'
#file_name_3 = '/home/zzchou/Documents/Codes/turner_distributions.txt'
#file_name_3 = '/home/zzchou/Documents/Codes/pierce_distributions.txt'
#file_name_4 = '/home/zzchou/Documents/Codes/epsztein_distributions.txt'
file_name_2 = '/home/zzchou/Documents/Codes/zmorph_tester_distributions.txt'
#file_name_3 = '/home/zzchou/Documents/Codes/gen2_zmorph_tester_distributions.txt'
#file_name_4 = '/home/zzchou/Documents/Codes/gen3_zmorph_tester_distributions.txt'
#file_name_2 = '/home/zzchou/Documents/Codes/gen4_zmorph_tester_distributions.txt'

legend_name_1 = file_name_1.split('/')[-1].split('_')[0]
legend_name_2 = file_name_2.split('/')[-1].split('_')[0]
#legend_name_3 = file_name_3.split('/')[-1].split('_')[0]
#legend_name_4 = file_name_4.split('/')[-1].split('_')[0]

distribution_names_1, distribution_type, distributions_1 = read_file(file_name_1)
distribution_names_2, distribution_type, distributions_2 = read_file(file_name_2)
#distribution_names_3, distribution_type, distributions_3 = read_file(file_name_3)
#distribution_names_4, distribution_type, distributions_4 = read_file(file_name_4)
basic_parameters = []
conditional_parameters = []

for x in range(0, len(distribution_type)):
    if distribution_type[x] == 'basic':
        basic_parameters.append(x)
    elif distribution_type[x] == 'conditional':
        conditional_parameters.append(x)
# plot_conditional_relationships(distribution_names_1, distributions_1, basic_parameters, conditional_parameters)
#compare_conditional_relationships(distribution_names_1, distributions_1, distributions_2, basic_parameters, conditional_parameters)

for i in range(0, min(len(distribution_names_1), len(distribution_names_2))):
    if distribution_names_1[i] == 'sholl bins':
        sholl_bins = distributions_1[i]
    elif distribution_names_1[i] == 'sholl counts':
        sholl_avgs_1 = [0]*len(sholl_bins)
        sholl_avgs_2 = [0]*len(sholl_bins)
#        sholl_avgs_3 = [0]*len(sholl_bins)
#        sholl_avgs_4 = [0]*len(sholl_bins)

        for j in range(0, len(sholl_bins)):
            for entry in distributions_1[i]:
                sholl_avgs_1[j]  = sholl_avgs_1[j] + entry[j]/len(distributions_1[i])
            for entry2 in distributions_2[i]:
                sholl_avgs_2[j] = sholl_avgs_2[j] + entry2[j]/len(distributions_2[i])
#            for entry3 in distributions_3[i]:
#                sholl_avgs_3[j] = sholl_avgs_3[j] + entry3[j]/len(distributions_3[i])
#            for entry4 in distributions_4[i]:
#                sholl_avgs_4[j] = sholl_avgs_4[j] + entry4[j]/len(distributions_4[i])

        sholl_avgs_1 = [x/max(sholl_avgs_1) for x in sholl_avgs_1]
        sholl_avgs_2 = [x/max(sholl_avgs_2) for x in sholl_avgs_2]
#        sholl_avgs_3 = [x/max(sholl_avgs_3) for x in sholl_avgs_3]
#        sholl_avgs_4 = [x/max(sholl_avgs_4) for x in sholl_avgs_4]
        
        print '\n\n\nSholl Analysis'
        plt.plot(sholl_bins, sholl_avgs_1, marker='o', label=legend_name_1)
        print sholl_bins, sholl_avgs_1
        plt.plot(sholl_bins, sholl_avgs_2, marker='o', label=legend_name_2)
#        plt.plot(sholl_bins, sholl_avgs_3, marker='o', label=legend_name_3)
#        plt.plot(sholl_bins, sholl_avgs_4, marker='o', label=legend_name_4)
        plt.xlabel('Distance from soma (um)', fontsize=18)
        plt.ylabel('Number of intersections', fontsize=18)
        plt.legend()
        plt.show()
    else:
        if i >=0 and i <= 200:
#            print '\n\n\n' + distribution_names_1[i]
            plot_2_distributions(distributions_1[i], distributions_2[i], distribution_names_1[i])
#            plot_3_distributions(distributions_1[i], distributions_2[i], distributions_3[i], distribution_names_1[i])
#            plot_4_distributions(distributions_1[i], distributions_2[i], distributions_3[i], distributions_4[i], distribution_names_1[i])

        
    