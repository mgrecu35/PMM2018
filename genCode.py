


def genCode(modName,tree, feature_names):
    left      = tree.tree_.children_left
    right     = tree.tree_.children_right
    threshold = tree.tree_.threshold
    features  = [feature_names[i] for i in tree.tree_.feature]
    value = tree.tree_.value
    s1=0
    scode= 'void %s('%modName
    for f in feature_names:
        scode+= 'float '+f+','
    scode+='float *value){\n'
    def recurse(left, right, threshold, features, node, s1,scode):
        if (threshold[node] != -2):
            scode+=\
                "if ( " + features[node] + " <= " + str(threshold[node]) + " ) {"
            if left[node] != -1:
                s1,scode=recurse (left, right, threshold, features,left[node],s1,scode)
                scode+= "} else {\n"
                if right[node] != -1:
                    s1,scode=recurse (left, right, threshold, features,right[node],s1,scode)
                    scode+= "}\n"
        else:
            scode+="*value="+str(value[node][0][0])+";\n"
            s1+=tree.tree_.n_node_samples[node]/(tree.tree_.n_node_samples[0]+0.
        return s1,scode
    s1,scode=recurse(left, right, threshold, features, 0, s1,scode)
    return scode+"}\n"
   

def genCodeEns(modName, feature_names,nfeatures):
    scode= 'void %s('%modName
    for f in feature_names:
        scode+= 'float '+f+','
    scode+='float *value){\n'
    scode+="int i;\n"
    scode+="float vali;\n"
    scode+="*value=0;\n"
    for i in range(nfeatures):
        fcall='%s%2.2i('%(modName,i)
        for f in feature_names:
            fcall+= f[1:]+','
        fcall+='&vali);\n'
        scode+=fcall
        scode+="*value=*value+vali;\n"
        #scode+="printf(\"%g\\n\",vali);"
    scode+="*value=*value/%i;\n"%nfeatures
    #scode+="printf(\"%g\\n\",*value);"
    scode+="}\n"
    return scode
