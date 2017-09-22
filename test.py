import regex

# the ?e flag grabs the best match (as few EDs as possible)
m=regex.findall("(?e)(GCTCCC){e<=1}", "CGCTCCA")  # means allow up to 1 error
print('Linker is ' + 'GCTCCC')
print(len('GCTCCC'))
print(m)
print(len(m[0]))

linker1 = 'TAGCCATCGCATTGC'
linker2 = 'TACCTCTAAGCTGAA'
regex_search = r".*([ACGT]{6})" + str(linker1) + r"([ACGT]{6,7})" \
                           + str(linker2) + r"([ACGT]{6})ACG([AGCT]{8}).?GACT"
match_obj = regex.match(regex_search, 'CTCTCAATTAGCCATCGCATTGCGGATAGGTACCTCTAAGCTGAAGGTGCTACGTGCCGCCCGACTTTT')
print('\n' + 'CTCTCAATTAGCCATCGCATTGCGGTAGGTACCTCTAAGCTGAAGGTGCTACGTGCCGCCCGACTTTT')
print(match_obj.group(2))
