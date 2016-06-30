
def seq_getter(filename1, outfile):
    """opens up a text file with gene names and distances to the nearest gene at the 5 prime
    and three prime direction. When given a list of wanted the scrpit will pull out the info
    for these genes of intesrest"""
    #text file with genes and distances - three columns
    f= open(filename1, "r")
    #outfile
    f_out= open(outfile, 'w')
    #assign the file contents to the variable data
    data = f.readlines()
    #load the data
    data1 = [line.rstrip("\n").split("\t") for line in (data)
             if line.strip() != "" and not line.startswith("g")]
    #to convert data 1 into a list of tuples.
    print data1[:10]
    location_data = [(str(s[0]), str(s[1]), str(s[2])) for s in (data1)]
    #genes of interest
    Dorsal = """GROS_g02635
GROS_g03475
GROS_g04620
GROS_g05241
GROS_g05352
GROS_g05354
GROS_g05707
GROS_g08876
GROS_g08879
GROS_g08949
GROS_g10807
GROS_g12028
GROS_g12196
GROS_g13195
GROS_g13646
GROS_g14123
GROS_g14125
GROS_g14157
GROS_g14158
GROS_g14179
GROS_g14180
GROS_g14194
GROS_g14228
GROS_g14234
GROS_g14236
GROS_g14275
GROS_g14306
    """.split()
    
    dpi_14 = """GROS_g08893
GROS_g12827
GROS_g00774
GROS_g01784
GROS_g02054
GROS_g03169
GROS_g04939
GROS_g06034
GROS_g06322
GROS_g07372
GROS_g07463
GROS_g08169
GROS_g08588
GROS_g08634
GROS_g08953
GROS_g08990
GROS_g09038
GROS_g09513
GROS_g09514
GROS_g09568
GROS_g10395
GROS_g10725
GROS_g10726
GROS_g11776
GROS_g11812
GROS_g12027
GROS_g12422
GROS_g12570
GROS_g12791
GROS_g13503
GROS_g13797
GROS_g14189
GROS_g14309
GROS_g06220
GROS_g08139
GROS_g14149
GROS_g02394
GROS_g02469
GROS_g02470
GROS_g05390
GROS_g05968
GROS_g06357
GROS_g08163
GROS_g09735
GROS_g10874
GROS_g10924
GROS_g11397
GROS_g11775
GROS_g11888
GROS_g12818
GROS_g13374
GROS_g14232
GROS_g03476
GROS_g06363
GROS_g06364
GROS_g09961
""".split()
    
    Dorsal_14 = """GROS_g00774
GROS_g01784
GROS_g02054
GROS_g03169
GROS_g04939
GROS_g06034
GROS_g06322
GROS_g07372
GROS_g07463
GROS_g08169
GROS_g08588
GROS_g08634
GROS_g08953
GROS_g08990
GROS_g09038
GROS_g09513
GROS_g09514
GROS_g09568
GROS_g10395
GROS_g10725
GROS_g10726
GROS_g11776
GROS_g11812
GROS_g12027
GROS_g12422
GROS_g12570
GROS_g12791
GROS_g13503
GROS_g13797
GROS_g14189
GROS_g14309
""".split()

    Subventral_14="""GROS_g03476
GROS_g06363
GROS_g06364
GROS_g09961""".split()

    J2 = """GROS_g02635
GROS_g03475
GROS_g04620
GROS_g05241
GROS_g05352
GROS_g05354
GROS_g05707
GROS_g08876
GROS_g08879
GROS_g08949
GROS_g10807
GROS_g12028
GROS_g12196
GROS_g13195
GROS_g13646
GROS_g14123
GROS_g14125
GROS_g14157
GROS_g14158
GROS_g14179
GROS_g14180
GROS_g14194
GROS_g14228
GROS_g14234
GROS_g14236
GROS_g14275
GROS_g14306
GROS_g02107
GROS_g03975
GROS_g00017
GROS_g04623
GROS_g05988
GROS_g06952
GROS_g07677
GROS_g08190
GROS_g09598
GROS_g14130
GROS_g14168
GROS_g14178
GROS_g14220
GROS_g14299
GROS_g05398
GROS_g05682
GROS_g06362
GROS_g07338
GROS_g07341
GROS_g11726
GROS_g08030
GROS_g08694
GROS_g12705
GROS_g04366
GROS_g04662
GROS_g04677
GROS_g05961
GROS_g05962
GROS_g05982
GROS_g06661
GROS_g07949
GROS_g07968
GROS_g10505
GROS_g10585
GROS_g11727
GROS_g12709
GROS_g14300""".split()

    Subventral_J2 = """GROS_g05398
GROS_g05682
GROS_g06362
GROS_g07338
GROS_g07341
GROS_g11726""".split()
    
    print >> f_out, 'geneid\tfiveprime\tthreeprime\tClass'
    name_set = set([])
    for i in location_data:
        gene_name = i[0].split(";")[0]
        if gene_name in Dorsal:
            if not gene_name in name_set:
                name_set.add(gene_name)
                fiveprime = i[1]
                threeprime = i[2]
                print >> f_out, '%s\t%s\t%s\tDorsal' %(gene_name, fiveprime,threeprime)
        if gene_name in dpi_14:
            if not gene_name in name_set:
                name_set.add(gene_name)
                fiveprime = i[1]
                threeprime = i[2]
                print >> f_out, '%s\t%s\t%s\tdpi_14' %(gene_name, fiveprime,threeprime)
        if gene_name in Dorsal_14:
            if not gene_name in name_set:
                name_set.add(gene_name)
                fiveprime = i[1]
                threeprime = i[2]
                print >> f_out, '%s\t%s\t%s\tDorsal_14' %(gene_name, fiveprime,threeprime)

        if gene_name in Subventral_14:
            if not gene_name in name_set:
                name_set.add(gene_name)
                fiveprime = i[1]
                threeprime = i[2]
                print >> f_out, '%s\t%s\t%s\tSubventral_14' %(gene_name, fiveprime,threeprime)

        if gene_name in J2:
            if not gene_name in name_set:
                name_set.add(gene_name)
                fiveprime = i[1]
                threeprime = i[2]
                print >> f_out, '%s\t%s\t%s\tJ2' %(gene_name, fiveprime,threeprime)

        if gene_name in Subventral_J2:
            if not gene_name in name_set:
                name_set.add(gene_name)
                fiveprime = i[1]
                threeprime = i[2]
                print >> f_out, '%s\t%s\t%s\tSubventral_J2' %(gene_name, fiveprime,threeprime)




            
        else:
            fiveprime = i[1]
            threeprime = i[2]
            print >> f_out, '%s\t%s\t%s\tOther' %(gene_name, fiveprime,threeprime)
    return True




seq_getter('distances_to_next_gene.txt', 'Gros_info_classes_for_heat_map001_20150429.txt')
print 'done'
