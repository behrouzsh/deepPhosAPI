import os
import json
import fnmatch
import sys
import urllib2
import pandas as pd
import numpy as np
import clustergrammer as clustergrammer
from flask import Flask, make_response, send_from_directory, request, render_template, url_for, redirect
from Model.models import db
from Service.ilincsService import IlincsSearch
from Service.make_clustergrammer import LoadFile
from scipy import stats





def get_valid_pyc_files():
    for import_path in sys.path:
        for root, dir_names, file_names in os.walk(import_path):
            for file_name in fnmatch.filter(file_names, '*.pyc'):
                pyc_file_path = os.path.join(root, file_name)
                py_file_path = pyc_file_path[:-1]

                if os.path.isfile(py_file_path):
                    yield pyc_file_path


def ensure_root():
    is_root = os.geteuid() == 0
    if not is_root:
        sys.exit('Must be run as root')


#def main():


def runserver():
    ensure_root()

    for pyc_file_path in get_valid_pyc_files():
        print 'removing', pyc_file_path
        os.remove(pyc_file_path)
    port = int(os.environ.get('PORT', 4000))
    app.run(host='0.0.0.0', port=port)

app = Flask(__name__)

# app.config['SQLALCHEMY_DATABASE_URI'] = 'postgresql://localhost/lincs2lincs'
# db.init_app(app)

app.secret_key = "development-key"

ilincsSearch = IlincsSearch()
loadFile = LoadFile()

# @app.route("/")
# @app.route("/home")
# @app.route("/about")
# @app.route("/clustergram")
# @app.route("/bootstrap")
# def basic_pages():
#     return make_response(open('static/index.html').read())




@app.route("/api/clust/P100/")
def make_clustergram_p100():
    petide2ptmAndGene = {
        "IYQY[+80]IQSR": "Q13627[Y+80@321](DYRK1A)",
        "TPKDS[+80]PGIPPSANAHQLFR": "P51812[S+80@369](RPS6KA3)",
        "RNS[+80]SEASSGDFLDLK": "Q9UK76[S+80@87](JPT1)",
        "LPLVPES[+80]PRR": "Q86WB0[S+80@321](ZC3HC1)",
        "ANAS[+80]PQKPLDLK": "Q9Y618[S+80@956](NCOR2)",
        "LENS[+80]PLGEALR": "Q9NX40[S+80@108](OCIAD1)",
        "ANS[+80]FVGTAQYVSPELLTEK": "O15530[S+80@241](PDPK1)",
        "TNPPTQKPPS[+80]PPMSGR": "Q8IZP0[S+80@183](ABI1)",
        "SNS[+80]LPHSAVSNAGSK": "Q8TBZ3[S+80@434](WDR20)",
        "VGS[+80]LDNVGHLPAGGAVK": "P27816[S+80@1073](MAP4)",
        "AAPEAS[+80]SPPASPLQHLLPGK": "Q96TA1[S+80@691](FAM129B)",
        "S[+122]DKPDM[+16]AEIEKFDK": "P62328[S+122@2][M+16@7](TMSB4X)",
        "S[+122]DKPDMAEIEKFDK": "P62328[S+122@2](TMSB4X)",
        "SLS[+80]LGDKEISR": "Q9UMZ2[S+80@1075](SYNRG)",
        "DLVQPDKPAS[+80]PK": "Q6PJT7[S+80@515](ZC3H14)",
        "SPS[+80]PAHLPDDPKVAEK": "Q92615[S+80@601](LARP4B)",
        "S[+80]IQDLTVTGTEPGQVSSR": "O43318[S+80@439](MAP3K7)",
        "IHS[+80]PIIR": "O60885[S+80@1117](BRD4)",
        "TFS[+80]LTEVR": "O95239[S+80@801](KIF4A)",
        "SLVGS[+80]WLK": "Q6ICG6[S+80@362](KIAA0930)",
        "S[+80]PPAPGLQPMR": "P15408[S+80@200](FOSL2)",
        "LAS[+80]PELER": "P17535[S+80@100](JUND)",
        "IGPLGLS[+80]PK": "P30050[S+80@38](RPL12)",
        "TPS[+80]IQPSLLPHAAPFAK": "P35658[S+80@1023](NUP214)",
        "HAS[+80]PILPITEFSDIPR": "P42167[S+80@306](TMPO)",
        "LIPGPLS[+80]PVAR": "P48634[S+80@1219](PRRC2A)",
        "LGM[+16]LS[+80]PEGTC[+57]K": "P49327[S+80@207][M+16@205][C+57@212](FASN)",
        "LGMLS[+80]PEGTC[+57]K": "P49327[S+80@207][C+57@212](FASN)",
        "ISNLS[+80]PEEEQGLWK": "Q5HYJ3[S+80@193](FAM76B)",
        "VSMPDVELNLKS[+80]PK": "Q09666[S+80@3426](AHNAK)",
        "S[+122]DNGELEDKPPAPPVR": "Q13177[S+122@2](PAK2)",
        "KAYS[+80]FC[+57]GTVEYM[+16]APEVVNR": "Q15418[S+80@221][M+16@229][C+57@223](RPS6KA1)",
        "KAYS[+80]FC[+57]GTVEYMAPEVVNR": "Q15418[S+80@221][C+57@223](RPS6KA1)",
        "NDS[+80]WGSFDLR": "Q7Z417[S+80@652](NUFIP2)",
        "LEVTEIVKPS[+80]PK": "Q7Z6E9[S+80@1179](RBBP6)",
        "YGS[+80]PPQRDPNWNGER": "O15234[S+80@265](CASC3)",
        "QDDS[+80]PPRPIIGPALPPGFIK": "Q8IXQ4[S+80@105](GPALPP1)",
        "SFS[+80]ADNFIGIQR": "Q8N7R7[S+80@344](CCNYL1)",
        "VLS[+80]PLIIK": "Q8NCN4[S+80@403](RNF169)",
        "AGS[+80]PDVLR": "Q8NDX6[S+80@44](ZNF740)",
        "LGPGRPLPTFPTSEC[+57]TS[+80]DVEPDTR": "Q8TDD1[S+80@75][C+57@73](DDX54)",
        "LAAPSVSHVS[+80]PR": "Q8WXE1[S+80@224](ATRIP)",
        "VDDDS[+80]LGEFPVTNSR": "Q92785[S+80@142](DPF2)",
        "NEEPVRS[+80]PERR": "Q92922[S+80@310](SMARCC1)",
        "LFIIRGS[+80]PQQIDHAK": "Q92945[S+80@480](KHSRP)",
        "S[+80]IEVENDFLPVEK": "Q96B97[S+80@230](SH3KBP1)",
        "TAPTLS[+80]PEHWK": "Q96JM3[S+80@405](CHAMP1)",
        "VLS[+80]PTAAKPSPFEGK": "Q96QC0[S+80@313](PPP1R10)",
        "SSDQPLTVPVS[+80]PK": "Q9ULW0[S+80@738](TPX2)",
        "FYETKEESYS[+80]PSKDR": "Q96T23[S+80@473](RSF1)",
        "SDS[+80]PENKYSDSTGHSK": "Q9BTA9[S+80@64](WAC)",
        "S[+80]IPLSIK": "Q9C0C9[S+80@515](UBE2O)",
        "RLS[+80]QSDEDVIR": "Q9H7D7[S+80@121](WDR26)",
        "ATS[+80]PVKSTTSITDAK": "Q9NQW6[S+80@295](ANLN)",
        "ALGS[+80]PTKQLLPC[+57]EMAC[+57]NEK": "Q9NR45[S+80@275][C+57@283][C+57@287](NANS)",
        "YLLGDAPVS[+80]PSSQK": "Q9NYB0[S+80@203](TERF2IP)",
        "ANS[+80]PEKPPEAGAAHKPR": "Q9UFC0[S+80@212](LRWD1)",
        "SEVQQPVHPKPLS[+80]PDSR": "Q9UHB6[S+80@362](LIMA1)",
        "ETPHS[+80]PGVEDAPIAK": "Q9UHB6[S+80@490](LIMA1)",
        "SQS[+80]PHYFR": "Q9UKJ3[S+80@1035](GPATCH8)",
        "DRS[+80]SPPPGYIPDELHQVAR": "Q9Y2U5[S+80@163](MAP3K2)",
        "SPALKS[+80]PLQSVVVR": "Q9Y2W1[S+80@253](THRAP3)",
        "AFGSGIDIKPGT[+80]PPIAGR": "Q9Y520[T+80@2673](PRRC2C)",
        "SFS[+80]SQRPVDR": "Q9Y520[S+80@1544](PRRC2C)",
        "VYT[+80]HEVVTLWYR": "P06493[T+80@161](CDK1)",
        "SST[+80]PLPTISSSAENTR": "P42167[T+80@160](TMPO)",
        "QIT[+80]MEELVR": "Q15149[T+80@4030](PLEC)",
        "TQLWASEPGT[+80]PPLPTSLPSQNPILK": "Q9BXP5[T+80@544](SRRT)",
        "ALPQT[+80]PRPR": "Q9UQ35[T+80@1492](SRRM2)",
        "SMS[+80]VDLSHIPLKDPLLFK": "A0JNW5[S+80@935](UHRF1BP1L)",
        "S[+80]PTGPSNSFLANMGGTVAHK": "Q96I25[S+80@222](RBM17)",
        "S[+80]LTAHSLLPLAEK": "Q86VI3[S+80@1424](IQGAP3)",
        "S[+80]FAGNLNTYKR": "Q01813[S+80@386](PFKP)",
        "HRPS[+80]PPATPPPK": "Q8IYB3[S+80@402](SRRM1)",
        "LHS[+80]APNLSDLHVVRPK": "O75385[S+80@556](ULK1)",
        "TLGRRDS[+80]SDDWEIPDGQITVGQR": "P15056[S+80@446](BRAF)",
        "A[+42]TTATM[+16]ATSGS[+80]AR": "P38919[S+80@12][M+16@7][A+42@2](EIF4A3)",
        "A[+42]TTATMATSGS[+80]AR": "P38919[S+80@12][A+42@2](EIF4A3)",
        "IHVSRS[+80]PTRPR": "Q499Z4[S+80@189](ZNF672)",
        "RPHS[+80]PEKAFSSNPVVR": "Q53F19[S+80@500](NCBP3)",
        "KPNIFYSGPAS[+80]PARPR": "Q6PL18[S+80@327](ATAD2)",
        "TEFLDLDNSPLSPPS[+80]PR": "Q8NCF5[S+80@204](NFATC2IP)",
        "QGSGRES[+80]PSLASR": "Q8WWM7[S+80@339](ATXN2L)",
        # "TQLWASEPGT[+80]PPLPTSLPSQNPILK": "Q8WWM7[S+80@339](SRRM2)",
        "LQS[+80]EPESIR": "P09496[S+80@105](CLTA)",
        "RLIS[+80]PYKK": "O14929[S+80@361](HAT1)",
        "LLEDS[+80]EESSEETVSR": "O60231[S+80@103](DHX16)",
        "S[+80]PPAPGLQPM[+16]R": "P15408[S+80@200][M+16@209](FOSL2)",
        "RRLS[+80]SLR": "P62753[S+80@235](RPS6)",
        "RLS[+80]ESQLSFRR": "Q96PK6[S+80@618](RBM14)",
        "RLS[+80]LPGLLSQVSPR": "Q96Q42[S+80@483](ALS2)",
        "SPDKPGGS[+80]PSASRR": "Q9Y3T9[S+80@56](NOC2L)",
        "HLPS[+80]PPTLDSIITEYLR": "Q9Y4B6[S+80@1000](DCAF1)",
        "ST[+80]FHAGQLR": "Q7KZI7[T+80@596](MARK2)",
        "S[+80]LTNSHLEKK": "Q9H2H9[S+80@52](SLC38A1)",
        "LQTPNT[+80]FPKR": "Q14978[T+80@610](NOLC1)",
        "QIT[+80]M[+16]EELVR": "Q15149[M+16@4031][T+80@4030](PLEC)",

        "T[+56]K[+56]QTAR": "P68431[me0K@5](HIST1H3A)",
        "T[+56]K[+70]QTAR": "P68431[meK@5](HIST1H3A)",
        "T[+56]K[+28]QTAR": "P68431[me2K@5](HIST1H3A)",
        "T[+56]K[+42]QTAR-me3K": "P68431[me3K@5](HIST1H3A)",
        "T[+56]K[+42]QTAR-aK": "P68431[aK@5](HIST1H3A)",
        "K[+112.1]STGGK[+56]APR": "P68431[me0K@10][a0K@15](HIST1H3A)",
        "K[+126.1]STGGK[+56]APR": "P68431[meK@10](HIST1H3A)",
        "K[+84.1]STGGK[+56]APR": "P68431[me2K@10](HIST1H3A)",
        "K[+98.1]STGGK[+56]APR": "P68431[me3K@10](HIST1H3A)",
        "K[+98]STGGK[+56]APR": "P68431[aK@10](HIST1H3A)",
        "K[+112.1]STGGK[+42]APR": "P68431[aK@15](HIST1H3A)",
        "K[+126.1]STGGK[+42]APR": "P68431[meK@10][aK@15](HIST1H3A)",
        "K[+84.1]STGGK[+42]APR": "P68431[me2K@10][aK@15](HIST1H3A)",
        "K[+98.1]STGGK[+42]APR": "P68431[me3K@10][aK@15](HIST1H3A)",
        "K[+98]STGGK[+42]APR": "P68431[aK@10][aK@15](HIST1H3A)",
        "K[+112.1]S[+80]TGGK[+56]APR": "P68431[pS@11](HIST1H3A)",
        "K[+126.1]S[+80]TGGK[+56]APR": "P68431[meK@10][pS@11](HIST1H3A)",
        "K[+84.1]S[+80]TGGK[+56]APR": "P68431[me2K@10][pS@11](HIST1H3A)",
        "K[+98.1]S[+80]TGGK[+56]APR": "P68431[me3K@10][pS@11](HIST1H3A)",
        "K[+98]S[+80]TGGK[+56]APR": "P68431[aK@10][pS@11](HIST1H3A)",
        "K[+112.1]S[+80]TGGK[+42]APR": "P68431[pS@11][aK@15](HIST1H3A)",
        "K[+126.1]S[+80]TGGK[+42]APR": "P68431[meK@10][pS@11][aK@15](HIST1H3A)",
        "K[+84.1]S[+80]TGGK[+42]APR": "P68431[me2K@10][pS@11][aK@15](HIST1H3A)",
        "K[+98.1]S[+80]TGGK[+42]APR": "P68431[me3K@10][pS@11][aK@15](HIST1H3A)",
        "K[+98]S[+80]TGGK[+42]APR": "P68431[aK@10][pS@11][aK@15](HIST1H3A)",

        "K[+112.1]QLATK[+56]AAR": "P68431[a0K@19][a0K@24](HIST1H3A)",
        "K[+98]QLATK[+56]AAR": "P68431[aK@19](HIST1H3A)",
        "K[+112.1]QLATK[+42]AAR": "P68431[aK@24](HIST1H3A)",
        "K[+98]QLATK[+42]AAR": "P68431[aK@19][aK@24](HIST1H3A)",

        "K[+226.1]QLATK[+56]AAR": "P68431[ubK@19][a0K@24](HIST1H3A)",
        "K[+112.1]QLATK[+170.1]AAR": "P68431[ubK@24](HIST1H3A)",
        "K[+112.1]SAPATGGVK[+56]K[+56]PHR": "P68431[me0K@28][me0K@37](HIST1H3A)",
        "K[+126.1]SAPATGGVK[+56]K[+56]PHR": "P68431[meK@28][me0K@37](HIST1H3A)",
        "K[+126.1]SAPATGGVK[+70]K[+56]PHR": "P68431[meK@28][meK@37](HIST1H3A)",
        "K[+126.1]SAPATGGVK[+28]K[+56]PHR": "P68431[meK@28][me2K@37](HIST1H3A)",
        "K[+126.1]SAPATGGVK[+42]K[+56]PHR": "P68431[meK@28][me3K@37](HIST1H3A)",
        "K[+84.1]SAPATGGVK[+56]K[+56]PHR": "P68431[me2K@28][me0K@37](HIST1H3A)",
        "K[+84.1]SAPATGGVK[+70]K[+56]PHR": "P68431[me2K@28][meK@37](HIST1H3A)",
        "K[+84.1]SAPATGGVK[+28]K[+56]PHR": "P68431[me2K@28][me2K@37](HIST1H3A)",
        "K[+84.1]SAPATGGVK[+42]K[+56]PHR": "P68431[me2K@28][me3K@37](HIST1H3A)",
        "K[+98.1]SAPATGGVK[+56]K[+56]PHR": "P68431[me3K@28][me0K@37](HIST1H3A)",
        "K[+98.1]SAPATGGVK[+70]K[+56]PHR": "P68431[me3K@28][meK@37](HIST1H3A)",
        "K[+98.1]SAPATGGVK[+28]K[+56]PHR": "P68431[me3K@28][me2K@37](HIST1H3A)",
        "K[+98.1]SAPATGGVK[+42]K[+56]PHR": "P68431[me3K@28][me3K@37](HIST1H3A)",
        "K[+98]SAPATGGVK[+56]K[+56]PHR": "P68431[aK@28][me0K@37](HIST1H3A)",
        "K[+98]SAPATGGVK[+70]K[+56]PHR": "P68431[aK@28][meK@37](HIST1H3A)",
        "K[+98]SAPATGGVK[+28]K[+56]PHR": "P68431[aK@28][me2K@37](HIST1H3A)",
        "K[+98]SAPATGGVK[+42]K[+56]PHR": "P68431[aK@28][me3K@37](HIST1H3A)",
        "K[+112.1]SAPSTGGVK[+56]K[+56]PHR": "P84243[me0K@28][me0K@37](H3F3A)",
        "Y[+56]RPGTVALR": "P68431-NORM(HIST1H3A)",
        "Y[+56]QK[+56]STELLIR": "P68431[me0K@57](HIST1H3A)",
        "E[+56]IAQDFK[+56]TDLR": "P68431[me0K@80](HIST1H3A)",
        "E[+56]IAQDFK[+70]TDLR": "P68431[meK@80](HIST1H3A)",
        "E[+56]IAQDFK[+28]TDLR": "P68431[me2K@80](HIST1H3A)",
        "K[+112.1]SAPATGGVK[+70]K[+56]PHR": "P68431[me0K@28][meK@37](HIST1H3A)",
        "K[+112.1]SAPATGGVK[+28]K[+56]PHR": "P68431[me0K@28][me2K@37](HIST1H3A)",
        "K[+112.1]SAPATGGVK[+42]K[+56]PHR": "P68431[me0K@28][me3K@37](HIST1H3A)",
        "Y[+56]QK[+42]STELLIR": "P68431[aK@57](HIST1H3A)",
        "G[+56]K[+56]GGK[+56]GLGK[+56]GGAK[+56]R": "P62805[a0K@6][a0K@13][a0K@17](HIST1H4A)",
        "G[+56]K[+42]GGK[+56]GLGK[+56]GGAK[+56]R": "P62805[aK@6](HIST1H4A)",
        "G[+56]K[+56]GGK[+56]GLGK[+42]GGAK[+56]R": "P62805[aK@13](HIST1H4A)",
        "G[+56]K[+56]GGK[+56]GLGK[+56]GGAK[+42]R": "P62805[aK@17](HIST1H4A)",
        "G[+56]K[+56]GGK[+42]GLGK[+42]GGAK[+56]R": "P62805[aK@9][aK@13](HIST1H4A)",
        "G[+56]K[+42]GGK[+42]GLGK[+56]GGAK[+56]R": "P62805[aK@6][aK@9](HIST1H4A)",
        "G[+56]K[+42]GGK[+56]GLGK[+56]GGAK[+42]R": "P62805[aK@6][aK@17](HIST1H4A)",
        "G[+56]K[+56]GGK[+56]GLGK[+42]GGAK[+42]R": "P62805[aK@13][aK@17](HIST1H4A)",
        "G[+56]K[+56]GGK[+42]GLGK[+42]GGAK[+42]R": "P62805[aK@9][aK@13][aK@17](HIST1H4A)",
        "G[+56]K[+42]GGK[+42]GLGK[+42]GGAK[+56]R": "P62805[aK@6][aK@9][aK@13](HIST1H4A)",
        "G[+56]K[+42]GGK[+42]GLGK[+56]GGAK[+42]R": "P62805[aK@6][aK@9][aK@17](HIST1H4A)",
        "G[+56]K[+42]GGK[+42]GLGK[+42]GGAK[+42]R": "P62805[aK@6][aK@9][aK@13][aK@17](HIST1H4A)",
        "G[+56]K[+56]GGK[+56]GLGK[+56]GGAK[+70]R": "P62805[meK@17](HIST1H4A)",
        "K[+112.1]VLR": "P62805[me0K@21](HIST1H4A)",
        "K[+126.1]VLR": "P62805[meK@21](HIST1H4A)",
        "K[+84.1]VLR": "P62805[me2K@21](HIST1H4A)",
        "K[+98.1]VLR": "P62805[me3K@21](HIST1H4A)",
        "D[+56]AVTYTEHAK[+56]R": "P62805-NORM(HIST1H4A)",
        "Y[+56]QK[+28]STELLIR": "P68431[me2K@57](HIST1H3A)"
    }
    P100_data = pd.read_csv('static/data/P100-all-plates-Level4.2018_08_07.processed.csv', sep=',', header=None)



    print P100_data.shape
    # print P100_data[0]
    print P100_data[0][5]
    print P100_data[12][5]
    col_cell_line = []
    col_pert_name = []
    col_pert_dose = []
    col_pert_time = []
    col_pert_type = []


    row_peptide = []
    row_protein = []
    row_gene = []
    row_shorthand = []

    for i in range(0, 96):
        row_peptide.append('Peptide: ' + str(P100_data[6][29 + i]))
        shorthand = petide2ptmAndGene[str(P100_data[6][29 + i])]
        row_shorthand.append('PTM Proteins: ' + shorthand)
        row_protein.append('Protein: ' + str(P100_data[11][29 + i]))
        row_gene.append('Gene: ' + str(P100_data[2][29 + i]))


    for column in range(21, P100_data.shape[1]):
        col_cell_line.append('Cell Line: ' + str(P100_data[column][1]))
        col_pert_dose.append('Dose: ' + str(P100_data[column][11]) + str(P100_data[column][12]))
        col_pert_time.append('Time: ' + str(P100_data[column][15]))
        col_pert_name.append('Perturnations: ' + str(P100_data[column][14]))
        if "_sg01" not in str(P100_data[column][14]) and "_sg02" not in str(P100_data[column][14]):
            col_pert_type.append('Type: Small Molecule')
        else:
            col_pert_type.append('Type: CRISPR')


    arrays_rows2 = [np.array(row_shorthand),np.array(row_peptide),np.array(row_protein),np.array(row_gene)]
    arrays_columns2 = [np.array(col_cell_line), np.array(col_pert_type), np.array(col_pert_name), np.array(col_pert_dose), np.array(col_pert_time)]


    tuples_rows2 = list(zip(*arrays_rows2))
    #print tuples_rows2
    tuples_columns2 = list(zip(*arrays_columns2))
    #print tuples_columns2
    rows_labels2 = pd.MultiIndex.from_tuples(tuples_rows2)
    columns_labels2 = pd.MultiIndex.from_tuples(tuples_columns2)


    a2 = np.matrix(P100_data.iloc[29:29+96, 21:P100_data.shape[1]].fillna(0))

    aa = a2.astype(np.float)
    print a2
    print a2.shape
    print columns_labels2.shape
    print rows_labels2.shape
    # {'nop': row1, 'o0p': row2, 'zaz': row3, 'zax': row4, 'oof': row5, 'oye': row6}
    # df = pd.DataFrame(a, index=rows_labels, columns=columns_labels)
    df2 = pd.DataFrame(aa, index=rows_labels2, columns=columns_labels2)
    #print df2
    P100_all_json = loadFile.make_json_from_txt(df2)
    json.dump(P100_all_json, open("static/data/P100_all_clustergram.json", 'w'))
    return P100_all_json



@app.route("/api/clust/GCP/aggregatedd/depricated/")
def make_clustergram_GCP_aggregated_depricated():
    petide2ptmAndGene = {
        "IYQY[+80]IQSR": "Q13627[Y+80@321](DYRK1A)",
        "TPKDS[+80]PGIPPSANAHQLFR": "P51812[S+80@369](RPS6KA3)",
        "RNS[+80]SEASSGDFLDLK": "Q9UK76[S+80@87](JPT1)",
        "LPLVPES[+80]PRR": "Q86WB0[S+80@321](ZC3HC1)",
        "ANAS[+80]PQKPLDLK": "Q9Y618[S+80@956](NCOR2)",
        "LENS[+80]PLGEALR": "Q9NX40[S+80@108](OCIAD1)",
        "ANS[+80]FVGTAQYVSPELLTEK": "O15530[S+80@241](PDPK1)",
        "TNPPTQKPPS[+80]PPMSGR": "Q8IZP0[S+80@183](ABI1)",
        "SNS[+80]LPHSAVSNAGSK": "Q8TBZ3[S+80@434](WDR20)",
        "VGS[+80]LDNVGHLPAGGAVK": "P27816[S+80@1073](MAP4)",
        "AAPEAS[+80]SPPASPLQHLLPGK": "Q96TA1[S+80@691](FAM129B)",
        "S[+122]DKPDM[+16]AEIEKFDK": "P62328[S+122@2][M+16@7](TMSB4X)",
        "S[+122]DKPDMAEIEKFDK": "P62328[S+122@2](TMSB4X)",
        "SLS[+80]LGDKEISR": "Q9UMZ2[S+80@1075](SYNRG)",
        "DLVQPDKPAS[+80]PK": "Q6PJT7[S+80@515](ZC3H14)",
        "SPS[+80]PAHLPDDPKVAEK": "Q92615[S+80@601](LARP4B)",
        "S[+80]IQDLTVTGTEPGQVSSR": "O43318[S+80@439](MAP3K7)",
        "IHS[+80]PIIR": "O60885[S+80@1117](BRD4)",
        "TFS[+80]LTEVR": "O95239[S+80@801](KIF4A)",
        "SLVGS[+80]WLK": "Q6ICG6[S+80@362](KIAA0930)",
        "S[+80]PPAPGLQPMR": "P15408[S+80@200](FOSL2)",
        "LAS[+80]PELER": "P17535[S+80@100](JUND)",
        "IGPLGLS[+80]PK": "P30050[S+80@38](RPL12)",
        "TPS[+80]IQPSLLPHAAPFAK": "P35658[S+80@1023](NUP214)",
        "HAS[+80]PILPITEFSDIPR": "P42167[S+80@306](TMPO)",
        "LIPGPLS[+80]PVAR": "P48634[S+80@1219](PRRC2A)",
        "LGM[+16]LS[+80]PEGTC[+57]K": "P49327[S+80@207][M+16@205][C+57@212](FASN)",
        "LGMLS[+80]PEGTC[+57]K": "P49327[S+80@207][C+57@212](FASN)",
        "ISNLS[+80]PEEEQGLWK": "Q5HYJ3[S+80@193](FAM76B)",
        "VSMPDVELNLKS[+80]PK": "Q09666[S+80@3426](AHNAK)",
        "S[+122]DNGELEDKPPAPPVR": "Q13177[S+122@2](PAK2)",
        "KAYS[+80]FC[+57]GTVEYM[+16]APEVVNR": "Q15418[S+80@221][M+16@229][C+57@223](RPS6KA1)",
        "KAYS[+80]FC[+57]GTVEYMAPEVVNR": "Q15418[S+80@221][C+57@223](RPS6KA1)",
        "NDS[+80]WGSFDLR": "Q7Z417[S+80@652](NUFIP2)",
        "LEVTEIVKPS[+80]PK": "Q7Z6E9[S+80@1179](RBBP6)",
        "YGS[+80]PPQRDPNWNGER": "O15234[S+80@265](CASC3)",
        "QDDS[+80]PPRPIIGPALPPGFIK": "Q8IXQ4[S+80@105](GPALPP1)",
        "SFS[+80]ADNFIGIQR": "Q8N7R7[S+80@344](CCNYL1)",
        "VLS[+80]PLIIK": "Q8NCN4[S+80@403](RNF169)",
        "AGS[+80]PDVLR": "Q8NDX6[S+80@44](ZNF740)",
        "LGPGRPLPTFPTSEC[+57]TS[+80]DVEPDTR": "Q8TDD1[S+80@75][C+57@73](DDX54)",
        "LAAPSVSHVS[+80]PR": "Q8WXE1[S+80@224](ATRIP)",
        "VDDDS[+80]LGEFPVTNSR": "Q92785[S+80@142](DPF2)",
        "NEEPVRS[+80]PERR": "Q92922[S+80@310](SMARCC1)",
        "LFIIRGS[+80]PQQIDHAK": "Q92945[S+80@480](KHSRP)",
        "S[+80]IEVENDFLPVEK": "Q96B97[S+80@230](SH3KBP1)",
        "TAPTLS[+80]PEHWK": "Q96JM3[S+80@405](CHAMP1)",
        "VLS[+80]PTAAKPSPFEGK": "Q96QC0[S+80@313](PPP1R10)",
        "SSDQPLTVPVS[+80]PK": "Q9ULW0[S+80@738](TPX2)",
        "FYETKEESYS[+80]PSKDR": "Q96T23[S+80@473](RSF1)",
        "SDS[+80]PENKYSDSTGHSK": "Q9BTA9[S+80@64](WAC)",
        "S[+80]IPLSIK": "Q9C0C9[S+80@515](UBE2O)",
        "RLS[+80]QSDEDVIR": "Q9H7D7[S+80@121](WDR26)",
        "ATS[+80]PVKSTTSITDAK": "Q9NQW6[S+80@295](ANLN)",
        "ALGS[+80]PTKQLLPC[+57]EMAC[+57]NEK": "Q9NR45[S+80@275][C+57@283][C+57@287](NANS)",
        "YLLGDAPVS[+80]PSSQK": "Q9NYB0[S+80@203](TERF2IP)",
        "ANS[+80]PEKPPEAGAAHKPR": "Q9UFC0[S+80@212](LRWD1)",
        "SEVQQPVHPKPLS[+80]PDSR": "Q9UHB6[S+80@362](LIMA1)",
        "ETPHS[+80]PGVEDAPIAK": "Q9UHB6[S+80@490](LIMA1)",
        "SQS[+80]PHYFR": "Q9UKJ3[S+80@1035](GPATCH8)",
        "DRS[+80]SPPPGYIPDELHQVAR": "Q9Y2U5[S+80@163](MAP3K2)",
        "SPALKS[+80]PLQSVVVR": "Q9Y2W1[S+80@253](THRAP3)",
        "AFGSGIDIKPGT[+80]PPIAGR": "Q9Y520[T+80@2673](PRRC2C)",
        "SFS[+80]SQRPVDR": "Q9Y520[S+80@1544](PRRC2C)",
        "VYT[+80]HEVVTLWYR": "P06493[T+80@161](CDK1)",
        "SST[+80]PLPTISSSAENTR": "P42167[T+80@160](TMPO)",
        "QIT[+80]MEELVR": "Q15149[T+80@4030](PLEC)",
        "TQLWASEPGT[+80]PPLPTSLPSQNPILK": "Q9BXP5[T+80@544](SRRT)",
        "ALPQT[+80]PRPR": "Q9UQ35[T+80@1492](SRRM2)",
        "SMS[+80]VDLSHIPLKDPLLFK": "A0JNW5[S+80@935](UHRF1BP1L)",
        "S[+80]PTGPSNSFLANMGGTVAHK": "Q96I25[S+80@222](RBM17)",
        "S[+80]LTAHSLLPLAEK": "Q86VI3[S+80@1424](IQGAP3)",
        "S[+80]FAGNLNTYKR": "Q01813[S+80@386](PFKP)",
        "HRPS[+80]PPATPPPK": "Q8IYB3[S+80@402](SRRM1)",
        "LHS[+80]APNLSDLHVVRPK": "O75385[S+80@556](ULK1)",
        "TLGRRDS[+80]SDDWEIPDGQITVGQR": "P15056[S+80@446](BRAF)",
        "A[+42]TTATM[+16]ATSGS[+80]AR": "P38919[S+80@12][M+16@7][A+42@2](EIF4A3)",
        "A[+42]TTATMATSGS[+80]AR": "P38919[S+80@12][A+42@2](EIF4A3)",
        "IHVSRS[+80]PTRPR": "Q499Z4[S+80@189](ZNF672)",
        "RPHS[+80]PEKAFSSNPVVR": "Q53F19[S+80@500](NCBP3)",
        "KPNIFYSGPAS[+80]PARPR": "Q6PL18[S+80@327](ATAD2)",
        "TEFLDLDNSPLSPPS[+80]PR": "Q8NCF5[S+80@204](NFATC2IP)",
        "QGSGRES[+80]PSLASR": "Q8WWM7[S+80@339](ATXN2L)",
        # "TQLWASEPGT[+80]PPLPTSLPSQNPILK": "Q8WWM7[S+80@339](SRRM2)",
        "LQS[+80]EPESIR": "P09496[S+80@105](CLTA)",
        "RLIS[+80]PYKK": "O14929[S+80@361](HAT1)",
        "LLEDS[+80]EESSEETVSR": "O60231[S+80@103](DHX16)",
        "S[+80]PPAPGLQPM[+16]R": "P15408[S+80@200][M+16@209](FOSL2)",
        "RRLS[+80]SLR": "P62753[S+80@235](RPS6)",
        "RLS[+80]ESQLSFRR": "Q96PK6[S+80@618](RBM14)",
        "RLS[+80]LPGLLSQVSPR": "Q96Q42[S+80@483](ALS2)",
        "SPDKPGGS[+80]PSASRR": "Q9Y3T9[S+80@56](NOC2L)",
        "HLPS[+80]PPTLDSIITEYLR": "Q9Y4B6[S+80@1000](DCAF1)",
        "ST[+80]FHAGQLR": "Q7KZI7[T+80@596](MARK2)",
        "S[+80]LTNSHLEKK": "Q9H2H9[S+80@52](SLC38A1)",
        "LQTPNT[+80]FPKR": "Q14978[T+80@610](NOLC1)",
        "QIT[+80]M[+16]EELVR": "Q15149[M+16@4031][T+80@4030](PLEC)",

        "T[+56]K[+56]QTAR": "P68431[me0K@5](HIST1H3A)",
        "T[+56]K[+70]QTAR": "P68431[meK@5](HIST1H3A)",
        "T[+56]K[+28]QTAR": "P68431[me2K@5](HIST1H3A)",
        "T[+56]K[+42]QTAR-me3K": "P68431[me3K@5](HIST1H3A)",
        "T[+56]K[+42]QTAR-aK": "P68431[aK@5](HIST1H3A)",
        "K[+112.1]STGGK[+56]APR": "P68431[me0K@10][a0K@15](HIST1H3A)",
        "K[+126.1]STGGK[+56]APR": "P68431[meK@10](HIST1H3A)",
        "K[+84.1]STGGK[+56]APR": "P68431[me2K@10](HIST1H3A)",
        "K[+98.1]STGGK[+56]APR": "P68431[me3K@10](HIST1H3A)",
        "K[+98]STGGK[+56]APR": "P68431[aK@10](HIST1H3A)",
        "K[+112.1]STGGK[+42]APR": "P68431[aK@15](HIST1H3A)",
        "K[+126.1]STGGK[+42]APR": "P68431[meK@10][aK@15](HIST1H3A)",
        "K[+84.1]STGGK[+42]APR": "P68431[me2K@10][aK@15](HIST1H3A)",
        "K[+98.1]STGGK[+42]APR": "P68431[me3K@10][aK@15](HIST1H3A)",
        "K[+98]STGGK[+42]APR": "P68431[aK@10][aK@15](HIST1H3A)",
        "K[+112.1]S[+80]TGGK[+56]APR": "P68431[pS@11](HIST1H3A)",
        "K[+126.1]S[+80]TGGK[+56]APR": "P68431[meK@10][pS@11](HIST1H3A)",
        "K[+84.1]S[+80]TGGK[+56]APR": "P68431[me2K@10][pS@11](HIST1H3A)",
        "K[+98.1]S[+80]TGGK[+56]APR": "P68431[me3K@10][pS@11](HIST1H3A)",
        "K[+98]S[+80]TGGK[+56]APR": "P68431[aK@10][pS@11](HIST1H3A)",
        "K[+112.1]S[+80]TGGK[+42]APR": "P68431[pS@11][aK@15](HIST1H3A)",
        "K[+126.1]S[+80]TGGK[+42]APR": "P68431[meK@10][pS@11][aK@15](HIST1H3A)",
        "K[+84.1]S[+80]TGGK[+42]APR": "P68431[me2K@10][pS@11][aK@15](HIST1H3A)",
        "K[+98.1]S[+80]TGGK[+42]APR": "P68431[me3K@10][pS@11][aK@15](HIST1H3A)",
        "K[+98]S[+80]TGGK[+42]APR": "P68431[aK@10][pS@11][aK@15](HIST1H3A)",

        "K[+112.1]QLATK[+56]AAR": "P68431[a0K@19][a0K@24](HIST1H3A)",
        "K[+98]QLATK[+56]AAR": "P68431[aK@19](HIST1H3A)",
        "K[+112.1]QLATK[+42]AAR": "P68431[aK@24](HIST1H3A)",
        "K[+98]QLATK[+42]AAR": "P68431[aK@19][aK@24](HIST1H3A)",

        "K[+226.1]QLATK[+56]AAR": "P68431[ubK@19][a0K@24](HIST1H3A)",
        "K[+112.1]QLATK[+170.1]AAR": "P68431[ubK@24](HIST1H3A)",
        "K[+112.1]SAPATGGVK[+56]K[+56]PHR": "P68431[me0K@28][me0K@37](HIST1H3A)",
        "K[+126.1]SAPATGGVK[+56]K[+56]PHR": "P68431[meK@28][me0K@37](HIST1H3A)",
        "K[+126.1]SAPATGGVK[+70]K[+56]PHR": "P68431[meK@28][meK@37](HIST1H3A)",
        "K[+126.1]SAPATGGVK[+28]K[+56]PHR": "P68431[meK@28][me2K@37](HIST1H3A)",
        "K[+126.1]SAPATGGVK[+42]K[+56]PHR": "P68431[meK@28][me3K@37](HIST1H3A)",
        "K[+84.1]SAPATGGVK[+56]K[+56]PHR": "P68431[me2K@28][me0K@37](HIST1H3A)",
        "K[+84.1]SAPATGGVK[+70]K[+56]PHR": "P68431[me2K@28][meK@37](HIST1H3A)",
        "K[+84.1]SAPATGGVK[+28]K[+56]PHR": "P68431[me2K@28][me2K@37](HIST1H3A)",
        "K[+84.1]SAPATGGVK[+42]K[+56]PHR": "P68431[me2K@28][me3K@37](HIST1H3A)",
        "K[+98.1]SAPATGGVK[+56]K[+56]PHR": "P68431[me3K@28][me0K@37](HIST1H3A)",
        "K[+98.1]SAPATGGVK[+70]K[+56]PHR": "P68431[me3K@28][meK@37](HIST1H3A)",
        "K[+98.1]SAPATGGVK[+28]K[+56]PHR": "P68431[me3K@28][me2K@37](HIST1H3A)",
        "K[+98.1]SAPATGGVK[+42]K[+56]PHR": "P68431[me3K@28][me3K@37](HIST1H3A)",
        "K[+98]SAPATGGVK[+56]K[+56]PHR": "P68431[aK@28][me0K@37](HIST1H3A)",
        "K[+98]SAPATGGVK[+70]K[+56]PHR": "P68431[aK@28][meK@37](HIST1H3A)",
        "K[+98]SAPATGGVK[+28]K[+56]PHR": "P68431[aK@28][me2K@37](HIST1H3A)",
        "K[+98]SAPATGGVK[+42]K[+56]PHR": "P68431[aK@28][me3K@37](HIST1H3A)",
        "K[+112.1]SAPSTGGVK[+56]K[+56]PHR": "P84243[me0K@28][me0K@37](H3F3A)",
        "Y[+56]RPGTVALR": "P68431-NORM(HIST1H3A)",
        "Y[+56]QK[+56]STELLIR": "P68431[me0K@57](HIST1H3A)",
        "E[+56]IAQDFK[+56]TDLR": "P68431[me0K@80](HIST1H3A)",
        "E[+56]IAQDFK[+70]TDLR": "P68431[meK@80](HIST1H3A)",
        "E[+56]IAQDFK[+28]TDLR": "P68431[me2K@80](HIST1H3A)",
        "K[+112.1]SAPATGGVK[+70]K[+56]PHR": "P68431[me0K@28][meK@37](HIST1H3A)",
        "K[+112.1]SAPATGGVK[+28]K[+56]PHR": "P68431[me0K@28][me2K@37](HIST1H3A)",
        "K[+112.1]SAPATGGVK[+42]K[+56]PHR": "P68431[me0K@28][me3K@37](HIST1H3A)",
        "Y[+56]QK[+42]STELLIR": "P68431[aK@57](HIST1H3A)",
        "G[+56]K[+56]GGK[+56]GLGK[+56]GGAK[+56]R": "P62805[a0K@6][a0K@13][a0K@17](HIST1H4A)",
        "G[+56]K[+42]GGK[+56]GLGK[+56]GGAK[+56]R": "P62805[aK@6](HIST1H4A)",
        "G[+56]K[+56]GGK[+56]GLGK[+42]GGAK[+56]R": "P62805[aK@13](HIST1H4A)",
        "G[+56]K[+56]GGK[+56]GLGK[+56]GGAK[+42]R": "P62805[aK@17](HIST1H4A)",
        "G[+56]K[+56]GGK[+42]GLGK[+42]GGAK[+56]R": "P62805[aK@9][aK@13](HIST1H4A)",
        "G[+56]K[+42]GGK[+42]GLGK[+56]GGAK[+56]R": "P62805[aK@6][aK@9](HIST1H4A)",
        "G[+56]K[+42]GGK[+56]GLGK[+56]GGAK[+42]R": "P62805[aK@6][aK@17](HIST1H4A)",
        "G[+56]K[+56]GGK[+56]GLGK[+42]GGAK[+42]R": "P62805[aK@13][aK@17](HIST1H4A)",
        "G[+56]K[+56]GGK[+42]GLGK[+42]GGAK[+42]R": "P62805[aK@9][aK@13][aK@17](HIST1H4A)",
        "G[+56]K[+42]GGK[+42]GLGK[+42]GGAK[+56]R": "P62805[aK@6][aK@9][aK@13](HIST1H4A)",
        "G[+56]K[+42]GGK[+42]GLGK[+56]GGAK[+42]R": "P62805[aK@6][aK@9][aK@17](HIST1H4A)",
        "G[+56]K[+42]GGK[+42]GLGK[+42]GGAK[+42]R": "P62805[aK@6][aK@9][aK@13][aK@17](HIST1H4A)",
        "G[+56]K[+56]GGK[+56]GLGK[+56]GGAK[+70]R": "P62805[meK@17](HIST1H4A)",
        "K[+112.1]VLR": "P62805[me0K@21](HIST1H4A)",
        "K[+126.1]VLR": "P62805[meK@21](HIST1H4A)",
        "K[+84.1]VLR": "P62805[me2K@21](HIST1H4A)",
        "K[+98.1]VLR": "P62805[me3K@21](HIST1H4A)",
        "D[+56]AVTYTEHAK[+56]R": "P62805-NORM(HIST1H4A)",
        "Y[+56]QK[+28]STELLIR": "P68431[me2K@57](HIST1H3A)"
    }
    GCP_data = pd.read_csv('static/data/GCP-all-plates-Level4.2018_08_07.processed.csv', sep=',', header=None)
    GCP_unique_complete = {}
    for column in range(25, GCP_data.shape[1]):
        # print (column)
        # identifier = str(GCP_data[column][10])  +"++"+str(GCP_data[column][3])  +"++"+ str(GCP_data[column][7]) +"++"+ str(GCP_data[column][11])
        GCP_identifier_complete = str(GCP_data[column][1]) + "++" + str(GCP_data[column][5]) + "++" + str(
            GCP_data[column][7]) + "++"  + str(GCP_data[column][8]) + "++" + str(GCP_data[column][9]) + "++" + \
                                   str(GCP_data[column][11]) + "++" + str(GCP_data[column][12]) + "++" + \
                                   str(GCP_data[column][13]) + "++" + str(GCP_data[column][14]) + "++" \
                                   + str(GCP_data[column][15]) + "++" + \
                                   str(GCP_data[column][16]) + "++" + str(GCP_data[column][17]) + "++" + \
                                   str(GCP_data[column][18]) + "++" + str(GCP_data[column][19]) + "++" + \
                                   str(GCP_data[column][20]) + "++" + str(GCP_data[column][21]) + "++" + \
                                   str(GCP_data[column][22]) + "++" + str(GCP_data[column][23]) + "++" + \
                                    str(GCP_data[column][24]) + "++" + str(GCP_data[column][25]) + "++" + \
                                    str(GCP_data[column][26]) + "++" + str(GCP_data[column][27])

        if GCP_identifier_complete not in GCP_unique_complete.keys():
            # print identifier_complete
            GCP_unique_complete[GCP_identifier_complete] = []
            GCP_unique_complete.get(GCP_identifier_complete).append(column)
        else:
            GCP_unique_complete.get(GCP_identifier_complete).append(column)


    iterator = 0
    for key in GCP_unique_complete:
        val = GCP_unique_complete.get(key)
        value = list(val)
        # print key
        for item in value:
            iterator += 1
        #     print item
        # print "------------------------"

    print iterator
    #print unique_complete_modified

    total = 0
    total_complete = 0


    print "size of unique_complete keys"
    print len(GCP_unique_complete.keys())
    print "size of unique_complete"
    print len(GCP_unique_complete)

    unique_sum = {}
    # unique_sum_modified = {}
    unique_pvalue = {}


    col_cell_line = []
    col_pert_name = []
    col_pert_dose = []
    col_pert_time = []
    col_pert_type = []

    my_data_con = []
    keyIter = 0
    first_loop = True
    for key in GCP_unique_complete:
        initial = [0.0]*79
        initial2 = [0.0] * 79
        initial_not_na = [0.0] * 79
        unique_sum[key] = initial
        unique_pvalue[key] = initial2
        val = GCP_unique_complete.get(key)
        key_splitted = key.split("++")
        print key_splitted
        col_cell_line.append('Cell Line: ' + key_splitted[0])
        col_pert_dose.append('Dose: ' + key_splitted[5] + key_splitted[6])
        col_pert_time.append('Time: ' + key_splitted[9] + key_splitted[10])
        col_pert_name.append('Perturnations: ' + key_splitted[8])
        if "_sg01" not in key_splitted[8] and "_sg02" not in key_splitted[8]:
            col_pert_type.append('Type: Small Molecule')
        else:
            col_pert_type.append('Type: CRISPR')

        item = list(val)
        # if len(item) == 1:
        #     print '111111   %d  =========== %s %s' % (keyIter, key, val)

        keyIter += 1

        value_list = 0
        for item_iter in range(0, 79):

            aux_list = []
            for col in item:
                aux_list.append(float(GCP_data[col][item_iter + 28]))
                if GCP_data[col][item_iter + 28] != "NA":
                    unique_sum.get(key)[item_iter] += float(GCP_data[col][item_iter + 28])
                    initial_not_na[item_iter] += 1
            if initial_not_na[item_iter] > 0:
                unique_sum.get(key)[item_iter] /= initial_not_na[item_iter]

            meanValue = np.mean(aux_list)
            average = sum(aux_list) / len(aux_list)
            # print "++++++++++++++++++"
            # print average
            # print meanValue
            # print unique_sum.get(key)[item_iter]
            if len(aux_list) > 1:

                variance2 = sum((average - value) ** 2 for value in aux_list) / (len(aux_list)-1)
        print "end of loop"

        x1 = np.array(unique_sum.get(key))
        print "end of loop1"
        x2 = x1[:, np.newaxis]
        print "end of loop2"
        my_data_final = np.asarray(x2)
        print "end of loop3"
        if first_loop:
            my_data_con = my_data_final
            first_loop = False

        else:
            my_data_con = np.hstack((my_data_con, my_data_final))
        print "end of loop4"
    print "======================"
    print "======================"
    print "======================"
    print my_data_con.shape




    row_peptide = []
    row_protein = []
    row_gene = []
    row_shorthand = []

    for i in range(0, 79):

        print str(GCP_data[4][28 + i])
        peptide = str(GCP_data[4][28 + i])
        if(peptide == 'T[+56]K[+42]QTAR'):
            print str(GCP_data[2][28 + i])
            if(str(GCP_data[2][28 + i]) == 'BI10006'):

                peptide = "T[+56]K[+42]QTAR-me3K"
            else:
                peptide = "T[+56]K[+42]QTAR-aK"
        print peptide
        row_peptide.append('Peptide: ' + peptide)

        shorthand = petide2ptmAndGene[peptide]

        print shorthand
        row_shorthand.append('PTM Proteins: ' + shorthand)

        row_protein.append('Protein: ' + str(GCP_data[9][28 + i]))
        row_gene.append('Gene: ' + str(GCP_data[6][28 + i]))


    arrays_rows2 = [np.array(row_shorthand),np.array(row_peptide),np.array(row_protein),np.array(row_gene)]
    arrays_columns2 = [np.array(col_pert_name), np.array(col_pert_time), np.array(col_pert_dose), np.array(col_pert_type),np.array(col_cell_line)]


    tuples_rows2 = list(zip(*arrays_rows2))
    #print tuples_rows2
    tuples_columns2 = list(zip(*arrays_columns2))
    #print tuples_columns2
    rows_labels2 = pd.MultiIndex.from_tuples(tuples_rows2)
    columns_labels2 = pd.MultiIndex.from_tuples(tuples_columns2)

    #my_data_con = my_data_con.fillna(0)
    print "after na to zero"
    a2 = np.nan_to_num(my_data_con)
    #a2 = my_data_con
    aa = a2.astype(np.float)
    print a2
    print a2.shape
    print columns_labels2.shape
    print rows_labels2.shape
    # {'nop': row1, 'o0p': row2, 'zaz': row3, 'zax': row4, 'oof': row5, 'oye': row6}
    #df = pd.DataFrame(a, index=rows_labels, columns=columns_labels)
    df2 = pd.DataFrame(aa, index=rows_labels2, columns=columns_labels2)
    #print df2
    GCP_all_json = loadFile.make_json_from_txt(df2)
    json.dump(GCP_all_json, open("static/data/GCP_aggregated_clustergram.json", 'w'))
    return GCP_all_json


@app.route("/api/test/<vars>")
def doSomething(vars):
    for var in vars.split(","):
        print var


@app.route("/api/clust/GCP/aggregated/")
def make_clustergram_GCP_aggregated():
    petide2ptmAndGene = {
        "IYQY[+80]IQSR": "Q13627[Y+80@321](DYRK1A)",
        "TPKDS[+80]PGIPPSANAHQLFR": "P51812[S+80@369](RPS6KA3)",
        "RNS[+80]SEASSGDFLDLK": "Q9UK76[S+80@87](JPT1)",
        "LPLVPES[+80]PRR": "Q86WB0[S+80@321](ZC3HC1)",
        "ANAS[+80]PQKPLDLK": "Q9Y618[S+80@956](NCOR2)",
        "LENS[+80]PLGEALR": "Q9NX40[S+80@108](OCIAD1)",
        "ANS[+80]FVGTAQYVSPELLTEK": "O15530[S+80@241](PDPK1)",
        "TNPPTQKPPS[+80]PPMSGR": "Q8IZP0[S+80@183](ABI1)",
        "SNS[+80]LPHSAVSNAGSK": "Q8TBZ3[S+80@434](WDR20)",
        "VGS[+80]LDNVGHLPAGGAVK": "P27816[S+80@1073](MAP4)",
        "AAPEAS[+80]SPPASPLQHLLPGK": "Q96TA1[S+80@691](FAM129B)",
        "S[+122]DKPDM[+16]AEIEKFDK": "P62328[S+122@2][M+16@7](TMSB4X)",
        "S[+122]DKPDMAEIEKFDK": "P62328[S+122@2](TMSB4X)",
        "SLS[+80]LGDKEISR": "Q9UMZ2[S+80@1075](SYNRG)",
        "DLVQPDKPAS[+80]PK": "Q6PJT7[S+80@515](ZC3H14)",
        "SPS[+80]PAHLPDDPKVAEK": "Q92615[S+80@601](LARP4B)",
        "S[+80]IQDLTVTGTEPGQVSSR": "O43318[S+80@439](MAP3K7)",
        "IHS[+80]PIIR": "O60885[S+80@1117](BRD4)",
        "TFS[+80]LTEVR": "O95239[S+80@801](KIF4A)",
        "SLVGS[+80]WLK": "Q6ICG6[S+80@362](KIAA0930)",
        "S[+80]PPAPGLQPMR": "P15408[S+80@200](FOSL2)",
        "LAS[+80]PELER": "P17535[S+80@100](JUND)",
        "IGPLGLS[+80]PK": "P30050[S+80@38](RPL12)",
        "TPS[+80]IQPSLLPHAAPFAK": "P35658[S+80@1023](NUP214)",
        "HAS[+80]PILPITEFSDIPR": "P42167[S+80@306](TMPO)",
        "LIPGPLS[+80]PVAR": "P48634[S+80@1219](PRRC2A)",
        "LGM[+16]LS[+80]PEGTC[+57]K": "P49327[S+80@207][M+16@205][C+57@212](FASN)",
        "LGMLS[+80]PEGTC[+57]K": "P49327[S+80@207][C+57@212](FASN)",
        "ISNLS[+80]PEEEQGLWK": "Q5HYJ3[S+80@193](FAM76B)",
        "VSMPDVELNLKS[+80]PK": "Q09666[S+80@3426](AHNAK)",
        "S[+122]DNGELEDKPPAPPVR": "Q13177[S+122@2](PAK2)",
        "KAYS[+80]FC[+57]GTVEYM[+16]APEVVNR": "Q15418[S+80@221][M+16@229][C+57@223](RPS6KA1)",
        "KAYS[+80]FC[+57]GTVEYMAPEVVNR": "Q15418[S+80@221][C+57@223](RPS6KA1)",
        "NDS[+80]WGSFDLR": "Q7Z417[S+80@652](NUFIP2)",
        "LEVTEIVKPS[+80]PK": "Q7Z6E9[S+80@1179](RBBP6)",
        "YGS[+80]PPQRDPNWNGER": "O15234[S+80@265](CASC3)",
        "QDDS[+80]PPRPIIGPALPPGFIK": "Q8IXQ4[S+80@105](GPALPP1)",
        "SFS[+80]ADNFIGIQR": "Q8N7R7[S+80@344](CCNYL1)",
        "VLS[+80]PLIIK": "Q8NCN4[S+80@403](RNF169)",
        "AGS[+80]PDVLR": "Q8NDX6[S+80@44](ZNF740)",
        "LGPGRPLPTFPTSEC[+57]TS[+80]DVEPDTR": "Q8TDD1[S+80@75][C+57@73](DDX54)",
        "LAAPSVSHVS[+80]PR": "Q8WXE1[S+80@224](ATRIP)",
        "VDDDS[+80]LGEFPVTNSR": "Q92785[S+80@142](DPF2)",
        "NEEPVRS[+80]PERR": "Q92922[S+80@310](SMARCC1)",
        "LFIIRGS[+80]PQQIDHAK": "Q92945[S+80@480](KHSRP)",
        "S[+80]IEVENDFLPVEK": "Q96B97[S+80@230](SH3KBP1)",
        "TAPTLS[+80]PEHWK": "Q96JM3[S+80@405](CHAMP1)",
        "VLS[+80]PTAAKPSPFEGK": "Q96QC0[S+80@313](PPP1R10)",
        "SSDQPLTVPVS[+80]PK": "Q9ULW0[S+80@738](TPX2)",
        "FYETKEESYS[+80]PSKDR": "Q96T23[S+80@473](RSF1)",
        "SDS[+80]PENKYSDSTGHSK": "Q9BTA9[S+80@64](WAC)",
        "S[+80]IPLSIK": "Q9C0C9[S+80@515](UBE2O)",
        "RLS[+80]QSDEDVIR": "Q9H7D7[S+80@121](WDR26)",
        "ATS[+80]PVKSTTSITDAK": "Q9NQW6[S+80@295](ANLN)",
        "ALGS[+80]PTKQLLPC[+57]EMAC[+57]NEK": "Q9NR45[S+80@275][C+57@283][C+57@287](NANS)",
        "YLLGDAPVS[+80]PSSQK": "Q9NYB0[S+80@203](TERF2IP)",
        "ANS[+80]PEKPPEAGAAHKPR": "Q9UFC0[S+80@212](LRWD1)",
        "SEVQQPVHPKPLS[+80]PDSR": "Q9UHB6[S+80@362](LIMA1)",
        "ETPHS[+80]PGVEDAPIAK": "Q9UHB6[S+80@490](LIMA1)",
        "SQS[+80]PHYFR": "Q9UKJ3[S+80@1035](GPATCH8)",
        "DRS[+80]SPPPGYIPDELHQVAR": "Q9Y2U5[S+80@163](MAP3K2)",
        "SPALKS[+80]PLQSVVVR": "Q9Y2W1[S+80@253](THRAP3)",
        "AFGSGIDIKPGT[+80]PPIAGR": "Q9Y520[T+80@2673](PRRC2C)",
        "SFS[+80]SQRPVDR": "Q9Y520[S+80@1544](PRRC2C)",
        "VYT[+80]HEVVTLWYR": "P06493[T+80@161](CDK1)",
        "SST[+80]PLPTISSSAENTR": "P42167[T+80@160](TMPO)",
        "QIT[+80]MEELVR": "Q15149[T+80@4030](PLEC)",
        "TQLWASEPGT[+80]PPLPTSLPSQNPILK": "Q9BXP5[T+80@544](SRRT)",
        "ALPQT[+80]PRPR": "Q9UQ35[T+80@1492](SRRM2)",
        "SMS[+80]VDLSHIPLKDPLLFK": "A0JNW5[S+80@935](UHRF1BP1L)",
        "S[+80]PTGPSNSFLANMGGTVAHK": "Q96I25[S+80@222](RBM17)",
        "S[+80]LTAHSLLPLAEK": "Q86VI3[S+80@1424](IQGAP3)",
        "S[+80]FAGNLNTYKR": "Q01813[S+80@386](PFKP)",
        "HRPS[+80]PPATPPPK": "Q8IYB3[S+80@402](SRRM1)",
        "LHS[+80]APNLSDLHVVRPK": "O75385[S+80@556](ULK1)",
        "TLGRRDS[+80]SDDWEIPDGQITVGQR": "P15056[S+80@446](BRAF)",
        "A[+42]TTATM[+16]ATSGS[+80]AR": "P38919[S+80@12][M+16@7][A+42@2](EIF4A3)",
        "A[+42]TTATMATSGS[+80]AR": "P38919[S+80@12][A+42@2](EIF4A3)",
        "IHVSRS[+80]PTRPR": "Q499Z4[S+80@189](ZNF672)",
        "RPHS[+80]PEKAFSSNPVVR": "Q53F19[S+80@500](NCBP3)",
        "KPNIFYSGPAS[+80]PARPR": "Q6PL18[S+80@327](ATAD2)",
        "TEFLDLDNSPLSPPS[+80]PR": "Q8NCF5[S+80@204](NFATC2IP)",
        "QGSGRES[+80]PSLASR": "Q8WWM7[S+80@339](ATXN2L)",
        # "TQLWASEPGT[+80]PPLPTSLPSQNPILK": "Q8WWM7[S+80@339](SRRM2)",
        "LQS[+80]EPESIR": "P09496[S+80@105](CLTA)",
        "RLIS[+80]PYKK": "O14929[S+80@361](HAT1)",
        "LLEDS[+80]EESSEETVSR": "O60231[S+80@103](DHX16)",
        "S[+80]PPAPGLQPM[+16]R": "P15408[S+80@200][M+16@209](FOSL2)",
        "RRLS[+80]SLR": "P62753[S+80@235](RPS6)",
        "RLS[+80]ESQLSFRR": "Q96PK6[S+80@618](RBM14)",
        "RLS[+80]LPGLLSQVSPR": "Q96Q42[S+80@483](ALS2)",
        "SPDKPGGS[+80]PSASRR": "Q9Y3T9[S+80@56](NOC2L)",
        "HLPS[+80]PPTLDSIITEYLR": "Q9Y4B6[S+80@1000](DCAF1)",
        "ST[+80]FHAGQLR": "Q7KZI7[T+80@596](MARK2)",
        "S[+80]LTNSHLEKK": "Q9H2H9[S+80@52](SLC38A1)",
        "LQTPNT[+80]FPKR": "Q14978[T+80@610](NOLC1)",
        "QIT[+80]M[+16]EELVR": "Q15149[M+16@4031][T+80@4030](PLEC)",

        "T[+56]K[+56]QTAR": "P68431[me0K@5](HIST1H3A)",
        "T[+56]K[+70]QTAR": "P68431[meK@5](HIST1H3A)",
        "T[+56]K[+28]QTAR": "P68431[me2K@5](HIST1H3A)",
        "T[+56]K[+42]QTAR-me3K": "P68431[me3K@5](HIST1H3A)",
        "T[+56]K[+42]QTAR-aK": "P68431[aK@5](HIST1H3A)",
        "K[+112.1]STGGK[+56]APR": "P68431[me0K@10][a0K@15](HIST1H3A)",
        "K[+126.1]STGGK[+56]APR": "P68431[meK@10](HIST1H3A)",
        "K[+84.1]STGGK[+56]APR": "P68431[me2K@10](HIST1H3A)",
        "K[+98.1]STGGK[+56]APR": "P68431[me3K@10](HIST1H3A)",
        "K[+98]STGGK[+56]APR": "P68431[aK@10](HIST1H3A)",
        "K[+112.1]STGGK[+42]APR": "P68431[aK@15](HIST1H3A)",
        "K[+126.1]STGGK[+42]APR": "P68431[meK@10][aK@15](HIST1H3A)",
        "K[+84.1]STGGK[+42]APR": "P68431[me2K@10][aK@15](HIST1H3A)",
        "K[+98.1]STGGK[+42]APR": "P68431[me3K@10][aK@15](HIST1H3A)",
        "K[+98]STGGK[+42]APR": "P68431[aK@10][aK@15](HIST1H3A)",
        "K[+112.1]S[+80]TGGK[+56]APR": "P68431[pS@11](HIST1H3A)",
        "K[+126.1]S[+80]TGGK[+56]APR": "P68431[meK@10][pS@11](HIST1H3A)",
        "K[+84.1]S[+80]TGGK[+56]APR": "P68431[me2K@10][pS@11](HIST1H3A)",
        "K[+98.1]S[+80]TGGK[+56]APR": "P68431[me3K@10][pS@11](HIST1H3A)",
        "K[+98]S[+80]TGGK[+56]APR": "P68431[aK@10][pS@11](HIST1H3A)",
        "K[+112.1]S[+80]TGGK[+42]APR": "P68431[pS@11][aK@15](HIST1H3A)",
        "K[+126.1]S[+80]TGGK[+42]APR": "P68431[meK@10][pS@11][aK@15](HIST1H3A)",
        "K[+84.1]S[+80]TGGK[+42]APR": "P68431[me2K@10][pS@11][aK@15](HIST1H3A)",
        "K[+98.1]S[+80]TGGK[+42]APR": "P68431[me3K@10][pS@11][aK@15](HIST1H3A)",
        "K[+98]S[+80]TGGK[+42]APR": "P68431[aK@10][pS@11][aK@15](HIST1H3A)",

        "K[+112.1]QLATK[+56]AAR": "P68431[a0K@19][a0K@24](HIST1H3A)",
        "K[+98]QLATK[+56]AAR": "P68431[aK@19](HIST1H3A)",
        "K[+112.1]QLATK[+42]AAR": "P68431[aK@24](HIST1H3A)",
        "K[+98]QLATK[+42]AAR": "P68431[aK@19][aK@24](HIST1H3A)",

        "K[+226.1]QLATK[+56]AAR": "P68431[ubK@19][a0K@24](HIST1H3A)",
        "K[+112.1]QLATK[+170.1]AAR": "P68431[ubK@24](HIST1H3A)",
        "K[+112.1]SAPATGGVK[+56]K[+56]PHR": "P68431[me0K@28][me0K@37](HIST1H3A)",
        "K[+126.1]SAPATGGVK[+56]K[+56]PHR": "P68431[meK@28][me0K@37](HIST1H3A)",
        "K[+126.1]SAPATGGVK[+70]K[+56]PHR": "P68431[meK@28][meK@37](HIST1H3A)",
        "K[+126.1]SAPATGGVK[+28]K[+56]PHR": "P68431[meK@28][me2K@37](HIST1H3A)",
        "K[+126.1]SAPATGGVK[+42]K[+56]PHR": "P68431[meK@28][me3K@37](HIST1H3A)",
        "K[+84.1]SAPATGGVK[+56]K[+56]PHR": "P68431[me2K@28][me0K@37](HIST1H3A)",
        "K[+84.1]SAPATGGVK[+70]K[+56]PHR": "P68431[me2K@28][meK@37](HIST1H3A)",
        "K[+84.1]SAPATGGVK[+28]K[+56]PHR": "P68431[me2K@28][me2K@37](HIST1H3A)",
        "K[+84.1]SAPATGGVK[+42]K[+56]PHR": "P68431[me2K@28][me3K@37](HIST1H3A)",
        "K[+98.1]SAPATGGVK[+56]K[+56]PHR": "P68431[me3K@28][me0K@37](HIST1H3A)",
        "K[+98.1]SAPATGGVK[+70]K[+56]PHR": "P68431[me3K@28][meK@37](HIST1H3A)",
        "K[+98.1]SAPATGGVK[+28]K[+56]PHR": "P68431[me3K@28][me2K@37](HIST1H3A)",
        "K[+98.1]SAPATGGVK[+42]K[+56]PHR": "P68431[me3K@28][me3K@37](HIST1H3A)",
        "K[+98]SAPATGGVK[+56]K[+56]PHR": "P68431[aK@28][me0K@37](HIST1H3A)",
        "K[+98]SAPATGGVK[+70]K[+56]PHR": "P68431[aK@28][meK@37](HIST1H3A)",
        "K[+98]SAPATGGVK[+28]K[+56]PHR": "P68431[aK@28][me2K@37](HIST1H3A)",
        "K[+98]SAPATGGVK[+42]K[+56]PHR": "P68431[aK@28][me3K@37](HIST1H3A)",
        "K[+112.1]SAPSTGGVK[+56]K[+56]PHR": "P84243[me0K@28][me0K@37](H3F3A)",
        "Y[+56]RPGTVALR": "P68431-NORM(HIST1H3A)",
        "Y[+56]QK[+56]STELLIR": "P68431[me0K@57](HIST1H3A)",
        "E[+56]IAQDFK[+56]TDLR": "P68431[me0K@80](HIST1H3A)",
        "E[+56]IAQDFK[+70]TDLR": "P68431[meK@80](HIST1H3A)",
        "E[+56]IAQDFK[+28]TDLR": "P68431[me2K@80](HIST1H3A)",
        "K[+112.1]SAPATGGVK[+70]K[+56]PHR": "P68431[me0K@28][meK@37](HIST1H3A)",
        "K[+112.1]SAPATGGVK[+28]K[+56]PHR": "P68431[me0K@28][me2K@37](HIST1H3A)",
        "K[+112.1]SAPATGGVK[+42]K[+56]PHR": "P68431[me0K@28][me3K@37](HIST1H3A)",
        "Y[+56]QK[+42]STELLIR": "P68431[aK@57](HIST1H3A)",
        "G[+56]K[+56]GGK[+56]GLGK[+56]GGAK[+56]R": "P62805[a0K@6][a0K@13][a0K@17](HIST1H4A)",
        "G[+56]K[+42]GGK[+56]GLGK[+56]GGAK[+56]R": "P62805[aK@6](HIST1H4A)",
        "G[+56]K[+56]GGK[+56]GLGK[+42]GGAK[+56]R": "P62805[aK@13](HIST1H4A)",
        "G[+56]K[+56]GGK[+56]GLGK[+56]GGAK[+42]R": "P62805[aK@17](HIST1H4A)",
        "G[+56]K[+56]GGK[+42]GLGK[+42]GGAK[+56]R": "P62805[aK@9][aK@13](HIST1H4A)",
        "G[+56]K[+42]GGK[+42]GLGK[+56]GGAK[+56]R": "P62805[aK@6][aK@9](HIST1H4A)",
        "G[+56]K[+42]GGK[+56]GLGK[+56]GGAK[+42]R": "P62805[aK@6][aK@17](HIST1H4A)",
        "G[+56]K[+56]GGK[+56]GLGK[+42]GGAK[+42]R": "P62805[aK@13][aK@17](HIST1H4A)",
        "G[+56]K[+56]GGK[+42]GLGK[+42]GGAK[+42]R": "P62805[aK@9][aK@13][aK@17](HIST1H4A)",
        "G[+56]K[+42]GGK[+42]GLGK[+42]GGAK[+56]R": "P62805[aK@6][aK@9][aK@13](HIST1H4A)",
        "G[+56]K[+42]GGK[+42]GLGK[+56]GGAK[+42]R": "P62805[aK@6][aK@9][aK@17](HIST1H4A)",
        "G[+56]K[+42]GGK[+42]GLGK[+42]GGAK[+42]R": "P62805[aK@6][aK@9][aK@13][aK@17](HIST1H4A)",
        "G[+56]K[+56]GGK[+56]GLGK[+56]GGAK[+70]R": "P62805[meK@17](HIST1H4A)",
        "K[+112.1]VLR": "P62805[me0K@21](HIST1H4A)",
        "K[+126.1]VLR": "P62805[meK@21](HIST1H4A)",
        "K[+84.1]VLR": "P62805[me2K@21](HIST1H4A)",
        "K[+98.1]VLR": "P62805[me3K@21](HIST1H4A)",
        "D[+56]AVTYTEHAK[+56]R": "P62805-NORM(HIST1H4A)",
        "Y[+56]QK[+28]STELLIR": "P68431[me2K@57](HIST1H3A)"
    }
    GCP_data = pd.read_csv('static/data/GCP-all-plates-Level4.2018_08_07.processed.csv', sep=',', header=None)
    GCP_unique_complete = {}
    print GCP_data.shape
    print "inside /api/clust/GCP/aggregated/"
    for column in range(25, GCP_data.shape[1]):
        print (column)
        # identifier = str(GCP_data[column][10])  +"++"+str(GCP_data[column][3])  +"++"+ str(GCP_data[column][7]) +"++"+ str(GCP_data[column][11])
        GCP_identifier_complete = str(GCP_data[column][1]) + "++" + str(GCP_data[column][5]) \
                                  + "++" + str(GCP_data[column][7]) + "++" + str(GCP_data[column][8]) + "++" + \
                                  str(GCP_data[column][11]) + "++" + str(GCP_data[column][12]) + "++" + \
                                   str(GCP_data[column][13]) + "++" + str(GCP_data[column][14]) + "++" \
                                   + str(GCP_data[column][15]) + "++" + \
                                   str(GCP_data[column][16]) + "++" + str(GCP_data[column][17]) + "++" + \
                                   str(GCP_data[column][18]) + "++" + str(GCP_data[column][19]) + "++" + \
                                   str(GCP_data[column][20]) + "++" + str(GCP_data[column][21]) + "++" + \
                                   str(GCP_data[column][22]) + "++" + str(GCP_data[column][23]) + "++" + \
                                   str(GCP_data[column][24]) + "++" + str(GCP_data[column][25]) + "++" + \
                                   str(GCP_data[column][26]) + "++" + str(GCP_data[column][27])

        if GCP_identifier_complete not in GCP_unique_complete.keys():
            # print identifier_complete
            GCP_unique_complete[GCP_identifier_complete] = []
            GCP_unique_complete.get(GCP_identifier_complete).append(column)
        else:
            GCP_unique_complete.get(GCP_identifier_complete).append(column)

    iterator = 0
    for key in GCP_unique_complete:
        val = GCP_unique_complete.get(key)
        value = list(val)
        # print key
        for item in value:
            iterator += 1
            #     print item
            # print "------------------------"

    print iterator
    # print unique_complete_modified

    total = 0
    total_complete = 0

    print "size of unique_complete keys"
    print len(GCP_unique_complete.keys())
    print "size of unique_complete"
    print len(GCP_unique_complete)

    unique_sum = {}
    # unique_sum_modified = {}
    unique_pvalue = {}

    col_cell_line = []
    col_pert_name = []
    col_pert_dose = []
    col_pert_time = []
    col_pert_type = []
    col_pert_signature = []

    my_data_con = []
    keyIter = 0
    first_loop = True

    # {"pertAnnotation": [{"LINCS_UNIQUE": "LINCSTP_0094", "pertVehicle": "DMSO", "pertName": "CPI-169", "cellId": "A375",
    #                    "pubchemCid": "71712226", "pertId": "BRD-A48881734", "pertType": "trt_cp", "pertTime": "3",
    #                   "pertDose": "1", "id": "P1", "lsmId": "LSM-42750"},
    #                 {"LINCS_UNIQUE": "LINCSTP_0236", "pertVehicle": "DMSO", "pertName": "CPI-169", "cellId": "MCF7",
    #                 "pubchemCid": "71712226", "pertId": "BRD-A48881734", "pertType": "trt_cp", "pertTime": "3",
    #                "pertDose": "1", "id": "P2", "lsmId": "LSM-42750"},
    #              {"LINCS_UNIQUE": "LINCSTP_0259", "pertVehicle": "DMSO", "pertName": "CPI-169", "cellId": "NPC",
    #              "pubchemCid": "71712226", "pertId": "BRD-A48881734", "pertType": "trt_cp", "pertTime": "3",
    #             "pertDose": "1", "id": "P3", "lsmId": "LSM-42750"},
    #           {"LINCS_UNIQUE": "LINCSTP_0116", "pertVehicle": "DMSO", "pertName": "CPI-169", "cellId": "A549",
    #           "pubchemCid": "71712226", "pertId": "BRD-A48881734", "pertType": "trt_cp", "pertTime": "3",
    #          "pertDose": "1", "id": "P4", "lsmId": "LSM-42750"},
    #        {"LINCS_UNIQUE": "LINCSTP_0275", "pertVehicle": "DMSO", "pertName": "CPI-169", "cellId": "PC3",
    #        "pubchemCid": "71712226", "pertId": "BRD-A48881734", "pertType": "trt_cp", "pertTime": "3",
    #       "pertDose": "1", "id": "P5", "lsmId": "LSM-42750"},
    #     {"LINCS_UNIQUE": "LINCSTP_0441", "pertVehicle": "DMSO", "pertName": "CPI-169", "cellId": "YAPC",
    #     "pubchemCid": "71712226", "pertId": "BRD-A48881734", "pertType": "trt_cp", "pertTime": "3",
    # "pertDose": "1", "id": "P6", "lsmId": "LSM-42750"}]}, "LINCSTP_0002":

    pertAnnotation = {}
    # this is for slicing the clustergrammer
    pertAnnotation_for_slicing = {}
    pertJson = {}
    unique_iterator = 10000
    keylist = GCP_unique_complete.keys()
    keylist.sort()
    for key in keylist:
        initial = [0.0] * 79
        initial2 = [0.0] * 79
        initial_not_na = [0.0] * 79
        signature_id = "LINCS_TP" + str(unique_iterator)
        unique_iterator += 1
        unique_sum[key] = initial
        unique_pvalue[key] = initial2
        val = GCP_unique_complete.get(key)
        key_splitted = key.split("++")
        print key_splitted
        col_cell_line.append('Cell Line: ' + key_splitted[0])
        col_pert_dose.append('Dose: ' + key_splitted[4] + key_splitted[5])
        col_pert_time.append('Time: ' + key_splitted[8] + key_splitted[9])
        col_pert_name.append(
            'Perturbations: ' + key_splitted[7] + "/" + key_splitted[4] + key_splitted[5] + "/" + key_splitted[8] +
            key_splitted[9])
        col_pert_signature.append('Signature ID: ' + signature_id)
        if "_sg01" not in key_splitted[7] and "_sg02" not in key_splitted[7]:
            cpType = 'Small Molecule'
        else:
            cpType = 'CRISPR'

        col_pert_type.append('Type: ' + cpType)
        compound_name = str(key_splitted[7]).upper()
        cp_metaData = {}
        cp_metaData_for_data_slicing = {}
        cp_metaData["Cell_line"] = key_splitted[0]
        cp_metaData["Dose"] = key_splitted[4] + " " + key_splitted[5]
        cp_metaData["Time"] = key_splitted[8] + " " + key_splitted[9]
        cp_metaData["Perturbations"] = compound_name
        cp_metaData["Type"] = cpType
        cp_metaData["Sig_id"] = signature_id

        cp_metaData_for_data_slicing["Cell_line"] = key_splitted[0]
        cp_metaData_for_data_slicing["Dose"] = key_splitted[4] + " " + key_splitted[5]
        cp_metaData_for_data_slicing["Time"] = key_splitted[8] + " " + key_splitted[9]
        cp_metaData_for_data_slicing["Perturbations"] = compound_name
        cp_metaData_for_data_slicing["Type"] = cpType
        cp_metaData_for_data_slicing["Sig_id"] = signature_id

        if compound_name not in pertAnnotation.keys():
            # print identifier_complete
            pertAnnotation[compound_name] = []
            pertAnnotation.get(compound_name).append(cp_metaData)
        else:
            pertAnnotation.get(compound_name).append(cp_metaData)

        item = list(val)
        # if len(item) == 1:
        #     print '111111   %d  =========== %s %s' % (keyIter, key, val)

        keyIter += 1

        value_list = 0
        sig2PeptideValue = {}
        for item_iter in range(0, 79):

            aux_list = []
            for col in item:

                if str(GCP_data[col][item_iter + 28]) != "nan":
                    aux_list.append(float(GCP_data[col][item_iter + 28]))
                    unique_sum.get(key)[item_iter] += float(GCP_data[col][item_iter + 28])
                    initial_not_na[item_iter] += 1
                else:
                    print cp_metaData
            if initial_not_na[item_iter] > 0:
                unique_sum.get(key)[item_iter] /= initial_not_na[item_iter]

            peptide = str(GCP_data[4][28 + item_iter])
            if (peptide == 'T[+56]K[+42]QTAR'):
                #print str(GCP_data[2][28 + item_iter])
                if (str(GCP_data[2][28 + item_iter]) == 'BI10006'):

                    peptide = "T[+56]K[+42]QTAR-me3K"
                else:
                    peptide = "T[+56]K[+42]QTAR-aK"
            shorthand = petide2ptmAndGene[peptide]
            sig2PeptideValue[shorthand] = unique_sum.get(key)[item_iter]

            meanValue = unique_sum.get(key)[item_iter]
            average = unique_sum.get(key)[item_iter]

            if len(aux_list) > 1:
                variance2 = sum((average - value) ** 2 for value in aux_list) / (len(aux_list) - 1)
                t_stat2 = meanValue / np.sqrt(variance2 / float(len(aux_list)))
                # p_val = stats.t.sf(abs(t_stat), len(item)-1)*2
                p_val2 = stats.t.sf(abs(t_stat2), len(aux_list) - 1) * 2
            else:
                p_val2 = 0.0
            unique_pvalue.get(key)[item_iter] = p_val2
            # variance = np.var(aux_list)

            # print aux_list
            # t_stat = meanValue/np.sqrt(variance/float(len(item)))













            if (average != unique_sum.get(key)[item_iter] or average != meanValue or meanValue != unique_sum.get(key)[
                item_iter]):
                print average
                print meanValue
                print unique_sum.get(key)[item_iter]
                print "++++++++++++++++++"
            if len(aux_list) > 1:
                variance2 = sum((average - value) ** 2 for value in aux_list) / (len(aux_list) - 1)
        # print "end of loop"
        cp_metaData_for_data_slicing["data"] = unique_sum.get(key)

        if compound_name not in pertAnnotation_for_slicing.keys():
            # print identifier_complete
            pertAnnotation_for_slicing[compound_name] = []
            pertAnnotation_for_slicing.get(compound_name).append(cp_metaData_for_data_slicing)
        else:
            pertAnnotation_for_slicing.get(compound_name).append(cp_metaData_for_data_slicing)

        x1 = np.array(unique_sum.get(key))
        # print "end of loop1"
        x2 = x1[:, np.newaxis]
        # print "end of loop2"
        my_data_final = np.asarray(x2)
        # print "end of loop3"
        if first_loop:
            my_data_con = my_data_final
            first_loop = False

        else:
            my_data_con = np.hstack((my_data_con, my_data_final))
        print "end of loop over unique"
        pertJson[signature_id] = sig2PeptideValue

    for keys in pertAnnotation:
        pertJson[keys] = pertAnnotation[keys]

    print "======================"
    print "======================"
    print "======================"
    print my_data_con.shape

    row_peptide = []
    row_protein = []
    row_gene = []
    row_shorthand = []

    for i in range(0, 79):

        peptide = str(GCP_data[4][28 + i])
        if(peptide == 'T[+56]K[+42]QTAR'):
            print str(GCP_data[2][28 + i])
            if(str(GCP_data[2][28 + i]) == 'BI10006'):

                peptide = "T[+56]K[+42]QTAR-me3K"
            else:
                peptide = "T[+56]K[+42]QTAR-aK"
        print peptide
        row_peptide.append('Peptide: ' + peptide)

        shorthand = petide2ptmAndGene[peptide]

        row_shorthand.append('PTM Proteins: ' + shorthand)
        row_protein.append('Protein: ' + str(GCP_data[9][28 + i]))
        row_gene.append('Gene: ' + str(GCP_data[6][28 + i]))


    arrays_rows2 = [np.array(row_shorthand), np.array(row_peptide), np.array(row_protein), np.array(row_gene)]
    arrays_columns2 = [np.array(col_pert_name), np.array(col_pert_signature), np.array(col_pert_time),
                       np.array(col_pert_dose), np.array(col_pert_type), np.array(col_cell_line)]
    np.save("static/data/GCP_rowMetaDat", arrays_rows2)
    tuples_rows2 = list(zip(*arrays_rows2))
    # print tuples_rows2
    tuples_columns2 = list(zip(*arrays_columns2))
    # print tuples_columns2
    rows_labels2 = pd.MultiIndex.from_tuples(tuples_rows2)
    columns_labels2 = pd.MultiIndex.from_tuples(tuples_columns2)

    # my_data_con = my_data_con.fillna(0)
    print "after na to zero"
    a2 = np.nan_to_num(my_data_con)
    # a2 = my_data_con
    aa = a2.astype(np.float)
    print a2
    print a2.shape
    print columns_labels2.shape
    print rows_labels2.shape
    # {'nop': row1, 'o0p': row2, 'zaz': row3, 'zax': row4, 'oof': row5, 'oye': row6}
    # df = pd.DataFrame(a, index=rows_labels, columns=columns_labels)
    df2 = pd.DataFrame(aa, index=rows_labels2, columns=columns_labels2)
    # print df2
    GCP_all_json = loadFile.make_json_from_txt(df2)
    json.dump(GCP_all_json, open("static/data/GCP_aggregated_clustergram.json", 'w'))
    json.dump(pertJson, open("static/data/GCP_processed_perturb_new.json", 'w'))
    json.dump(pertAnnotation_for_slicing, open("static/data/GCP_processed_perturb_for_clustergramm_slicing.json", 'w'))

    # str(GCP_data[column][1]) + "++" + str(GCP_data[column][3]) + "++" + str(
    #     GCP_data[column][4]) + "++" + str(GCP_data[column][7]) + "++" + \
    # str(GCP_data[column][8]) + "++" + \
    # str(GCP_data[column][11]) + "++" + str(GCP_data[column][12]) + "++" + \
    # str(GCP_data[column][13]) + "++" + str(GCP_data[column][14]) + "++" \
    # + str(GCP_data[column][15]) + "++" + \
    # str(GCP_data[column][16]) + "++" + str(GCP_data[column][17]) + "++" + \
    # str(GCP_data[column][18]) + "++" + str(GCP_data[column][19]) + "++" + \
    # str(GCP_data[column][20]) + "++" + str(GCP_data[column][21]) + "++" + \
    # str(GCP_data[column][22]) + "++" + str(GCP_data[column][23]) + "++" + \
    # str(GCP_data[column][24]) + "++" + str(GCP_data[column][25]) + "++" + \
    # str(GCP_data[column][26]) + "++" + str(GCP_data[column][27]) + "++" + \
    # str(GCP_data[column][28])

    # making p-value file
    first_col = []
    first_col.append("")
    first_col.append("LINCS_Signature_ID")
    first_col.append(str(GCP_data[0][1]))
    first_col.append(str(GCP_data[0][5]))
    first_col.append(str(GCP_data[0][7]))
    first_col.append(str(GCP_data[0][8]))

    first_col.append(str(GCP_data[0][11]))
    first_col.append(str(GCP_data[0][12]))
    first_col.append(str(GCP_data[0][13]))
    first_col.append(str(GCP_data[0][14]))
    first_col.append(str(GCP_data[0][15]))
    first_col.append(str(GCP_data[0][16]))
    first_col.append(str(GCP_data[0][17]))
    first_col.append(str(GCP_data[0][18]))
    first_col.append(str(GCP_data[0][19]))
    first_col.append(str(GCP_data[0][20]))
    first_col.append(str(GCP_data[0][21]))
    first_col.append(str(GCP_data[0][22]))
    first_col.append(str(GCP_data[0][23]))
    first_col.append(str(GCP_data[0][24]))
    first_col.append(str(GCP_data[0][25]))
    first_col.append(str(GCP_data[0][26]))
    first_col.append(str(GCP_data[0][27]))


    for i in range(0, 79):
        first_col.append(str(GCP_data[0][28 + i]))

    x1 = np.array(first_col)
    x2 = x1[:, np.newaxis]
    my_data_final = np.asarray(x2)

    for col in range(1, 11):
        ith_col = []
        ith_col.append(str(GCP_data[col][0]))
        ith_col.append(str(""))
        ith_col.append(str(""))
        ith_col.append(str(""))
        ith_col.append(str(""))
        ith_col.append(str(""))
        ith_col.append(str(""))
        ith_col.append(str(""))
        ith_col.append(str(""))
        ith_col.append(str(""))
        ith_col.append(str(""))
        ith_col.append(str(""))
        ith_col.append(str(""))
        ith_col.append(str(""))
        ith_col.append(str(""))
        ith_col.append(str(""))
        ith_col.append(str(""))
        ith_col.append(str(""))
        ith_col.append(str(""))
        ith_col.append(str(""))
        ith_col.append(str(""))
        ith_col.append(str(""))
        ith_col.append(str(""))



        for i in range(0, 79):
            ith_col.append(str(GCP_data[col][28 + i]))

        x1 = np.array(ith_col)
        x2 = x1[:, np.newaxis]
        col_data_final = np.asarray(x2)
        my_data_final = np.hstack((my_data_final, col_data_final))

    # Adding PTMProtein
    ith_col = []
    ith_col.append(str("ShortHand"))
    ith_col.append(str(""))
    ith_col.append(str(""))
    ith_col.append(str(""))
    ith_col.append(str(""))
    ith_col.append(str(""))
    ith_col.append(str(""))
    ith_col.append(str(""))
    ith_col.append(str(""))
    ith_col.append(str(""))
    ith_col.append(str(""))
    ith_col.append(str(""))
    ith_col.append(str(""))
    ith_col.append(str(""))
    ith_col.append(str(""))
    ith_col.append(str(""))
    ith_col.append(str(""))
    ith_col.append(str(""))
    ith_col.append(str(""))
    ith_col.append(str(""))
    ith_col.append(str(""))
    ith_col.append(str(""))
    ith_col.append(str(""))



    for i in range(0, 79):
        peptide = str(GCP_data[4][28 + i])
        if(peptide == 'T[+56]K[+42]QTAR'):
            print str(GCP_data[2][28 + i])
            if(str(GCP_data[2][28 + i]) == 'BI10006'):

                peptide = "T[+56]K[+42]QTAR-me3K"
            else:
                peptide = "T[+56]K[+42]QTAR-aK"


        shorthand = petide2ptmAndGene[peptide]

        ith_col.append(shorthand)
    x1 = np.array(ith_col)
    x2 = x1[:, np.newaxis]
    col_data_final = np.asarray(x2)
    my_data_final = np.hstack((my_data_final, col_data_final))


    my_data_final2 = my_data_final
    my_data_final3 = my_data_final
    my_pvalue_con = my_data_final2
    my_average_con = my_data_final3
    print "my_data_final.shape"
    print my_data_final.shape
    unique_iterator = 10000
    keylist = GCP_unique_complete.keys()
    keylist.sort()
    for key in keylist:

        signature_id = "LINCS_TP" + str(unique_iterator)
        unique_iterator += 1

        key_splitted = key.split("++")

        key_list1 = []
        key_list1.append("")
        key_list1.append(signature_id)

        for list_member in key_splitted:
            key_list1.append(list_member)

        # change list to array
        key_arr1 = np.array(key_list1)

        # This is for making data n by 1 instead of 1 by n
        key_arr2 = key_arr1[:, np.newaxis]

        x1 = np.array(unique_sum.get(key))
        x2 = x1[:, np.newaxis]

        my_data = np.asarray(x2)
        col_data = np.vstack((key_arr2, my_data))

        # Repeating for p-values
        x3 = np.array(unique_pvalue.get(key))
        x4 = x3[:, np.newaxis]

        my_data_pvalue = np.asarray(x4)
        col_data_pvalue = np.vstack((key_arr2, my_data_pvalue))

        my_average_con = np.hstack((my_average_con, col_data))
        my_pvalue_con = np.hstack((my_pvalue_con, col_data_pvalue))

    print "my_data.shape"
    print my_pvalue_con.shape

    np.savetxt("static/data/GCP_average.txt", my_average_con, delimiter='\t', fmt="%s")
    np.savetxt("static/data/GCP_pValue.txt", my_pvalue_con, delimiter='\t', fmt="%s")

    return GCP_all_json

@app.route("/api/clust/P100/aggregated/")
def make_clustergram_p100_aggregated():
    petide2ptmAndGene = {
        "IYQY[+80]IQSR": "Q13627[Y+80@321](DYRK1A)",
        "TPKDS[+80]PGIPPSANAHQLFR": "P51812[S+80@369](RPS6KA3)",
        "RNS[+80]SEASSGDFLDLK": "Q9UK76[S+80@87](JPT1)",
        "LPLVPES[+80]PRR": "Q86WB0[S+80@321](ZC3HC1)",
        "ANAS[+80]PQKPLDLK": "Q9Y618[S+80@956](NCOR2)",
        "LENS[+80]PLGEALR": "Q9NX40[S+80@108](OCIAD1)",
        "ANS[+80]FVGTAQYVSPELLTEK": "O15530[S+80@241](PDPK1)",
        "TNPPTQKPPS[+80]PPMSGR": "Q8IZP0[S+80@183](ABI1)",
        "SNS[+80]LPHSAVSNAGSK": "Q8TBZ3[S+80@434](WDR20)",
        "VGS[+80]LDNVGHLPAGGAVK": "P27816[S+80@1073](MAP4)",
        "AAPEAS[+80]SPPASPLQHLLPGK": "Q96TA1[S+80@691](FAM129B)",
        "S[+122]DKPDM[+16]AEIEKFDK": "P62328[S+122@2][M+16@7](TMSB4X)",
        "S[+122]DKPDMAEIEKFDK": "P62328[S+122@2](TMSB4X)",
        "SLS[+80]LGDKEISR": "Q9UMZ2[S+80@1075](SYNRG)",
        "DLVQPDKPAS[+80]PK": "Q6PJT7[S+80@515](ZC3H14)",
        "SPS[+80]PAHLPDDPKVAEK": "Q92615[S+80@601](LARP4B)",
        "S[+80]IQDLTVTGTEPGQVSSR": "O43318[S+80@439](MAP3K7)",
        "IHS[+80]PIIR": "O60885[S+80@1117](BRD4)",
        "TFS[+80]LTEVR": "O95239[S+80@801](KIF4A)",
        "SLVGS[+80]WLK": "Q6ICG6[S+80@362](KIAA0930)",
        "S[+80]PPAPGLQPMR": "P15408[S+80@200](FOSL2)",
        "LAS[+80]PELER": "P17535[S+80@100](JUND)",
        "IGPLGLS[+80]PK": "P30050[S+80@38](RPL12)",
        "TPS[+80]IQPSLLPHAAPFAK": "P35658[S+80@1023](NUP214)",
        "HAS[+80]PILPITEFSDIPR": "P42167[S+80@306](TMPO)",
        "LIPGPLS[+80]PVAR": "P48634[S+80@1219](PRRC2A)",
        "LGM[+16]LS[+80]PEGTC[+57]K": "P49327[S+80@207][M+16@205][C+57@212](FASN)",
        "LGMLS[+80]PEGTC[+57]K": "P49327[S+80@207][C+57@212](FASN)",
        "ISNLS[+80]PEEEQGLWK": "Q5HYJ3[S+80@193](FAM76B)",
        "VSMPDVELNLKS[+80]PK": "Q09666[S+80@3426](AHNAK)",
        "S[+122]DNGELEDKPPAPPVR": "Q13177[S+122@2](PAK2)",
        "KAYS[+80]FC[+57]GTVEYM[+16]APEVVNR": "Q15418[S+80@221][M+16@229][C+57@223](RPS6KA1)",
        "KAYS[+80]FC[+57]GTVEYMAPEVVNR": "Q15418[S+80@221][C+57@223](RPS6KA1)",
        "NDS[+80]WGSFDLR": "Q7Z417[S+80@652](NUFIP2)",
        "LEVTEIVKPS[+80]PK": "Q7Z6E9[S+80@1179](RBBP6)",
        "YGS[+80]PPQRDPNWNGER": "O15234[S+80@265](CASC3)",
        "QDDS[+80]PPRPIIGPALPPGFIK": "Q8IXQ4[S+80@105](GPALPP1)",
        "SFS[+80]ADNFIGIQR": "Q8N7R7[S+80@344](CCNYL1)",
        "VLS[+80]PLIIK": "Q8NCN4[S+80@403](RNF169)",
        "AGS[+80]PDVLR": "Q8NDX6[S+80@44](ZNF740)",
        "LGPGRPLPTFPTSEC[+57]TS[+80]DVEPDTR": "Q8TDD1[S+80@75][C+57@73](DDX54)",
        "LAAPSVSHVS[+80]PR": "Q8WXE1[S+80@224](ATRIP)",
        "VDDDS[+80]LGEFPVTNSR": "Q92785[S+80@142](DPF2)",
        "NEEPVRS[+80]PERR": "Q92922[S+80@310](SMARCC1)",
        "LFIIRGS[+80]PQQIDHAK": "Q92945[S+80@480](KHSRP)",
        "S[+80]IEVENDFLPVEK": "Q96B97[S+80@230](SH3KBP1)",
        "TAPTLS[+80]PEHWK": "Q96JM3[S+80@405](CHAMP1)",
        "VLS[+80]PTAAKPSPFEGK": "Q96QC0[S+80@313](PPP1R10)",
        "SSDQPLTVPVS[+80]PK": "Q9ULW0[S+80@738](TPX2)",
        "FYETKEESYS[+80]PSKDR": "Q96T23[S+80@473](RSF1)",
        "SDS[+80]PENKYSDSTGHSK": "Q9BTA9[S+80@64](WAC)",
        "S[+80]IPLSIK": "Q9C0C9[S+80@515](UBE2O)",
        "RLS[+80]QSDEDVIR": "Q9H7D7[S+80@121](WDR26)",
        "ATS[+80]PVKSTTSITDAK": "Q9NQW6[S+80@295](ANLN)",
        "ALGS[+80]PTKQLLPC[+57]EMAC[+57]NEK": "Q9NR45[S+80@275][C+57@283][C+57@287](NANS)",
        "YLLGDAPVS[+80]PSSQK": "Q9NYB0[S+80@203](TERF2IP)",
        "ANS[+80]PEKPPEAGAAHKPR": "Q9UFC0[S+80@212](LRWD1)",
        "SEVQQPVHPKPLS[+80]PDSR": "Q9UHB6[S+80@362](LIMA1)",
        "ETPHS[+80]PGVEDAPIAK": "Q9UHB6[S+80@490](LIMA1)",
        "SQS[+80]PHYFR": "Q9UKJ3[S+80@1035](GPATCH8)",
        "DRS[+80]SPPPGYIPDELHQVAR": "Q9Y2U5[S+80@163](MAP3K2)",
        "SPALKS[+80]PLQSVVVR": "Q9Y2W1[S+80@253](THRAP3)",
        "AFGSGIDIKPGT[+80]PPIAGR": "Q9Y520[T+80@2673](PRRC2C)",
        "SFS[+80]SQRPVDR": "Q9Y520[S+80@1544](PRRC2C)",
        "VYT[+80]HEVVTLWYR": "P06493[T+80@161](CDK1)",
        "SST[+80]PLPTISSSAENTR": "P42167[T+80@160](TMPO)",
        "QIT[+80]MEELVR": "Q15149[T+80@4030](PLEC)",
        "TQLWASEPGT[+80]PPLPTSLPSQNPILK": "Q9BXP5[T+80@544](SRRT)",
        "ALPQT[+80]PRPR": "Q9UQ35[T+80@1492](SRRM2)",
        "SMS[+80]VDLSHIPLKDPLLFK": "A0JNW5[S+80@935](UHRF1BP1L)",
        "S[+80]PTGPSNSFLANMGGTVAHK": "Q96I25[S+80@222](RBM17)",
        "S[+80]LTAHSLLPLAEK": "Q86VI3[S+80@1424](IQGAP3)",
        "S[+80]FAGNLNTYKR": "Q01813[S+80@386](PFKP)",
        "HRPS[+80]PPATPPPK": "Q8IYB3[S+80@402](SRRM1)",
        "LHS[+80]APNLSDLHVVRPK": "O75385[S+80@556](ULK1)",
        "TLGRRDS[+80]SDDWEIPDGQITVGQR": "P15056[S+80@446](BRAF)",
        "A[+42]TTATM[+16]ATSGS[+80]AR": "P38919[S+80@12][M+16@7][A+42@2](EIF4A3)",
        "A[+42]TTATMATSGS[+80]AR": "P38919[S+80@12][A+42@2](EIF4A3)",
        "IHVSRS[+80]PTRPR": "Q499Z4[S+80@189](ZNF672)",
        "RPHS[+80]PEKAFSSNPVVR": "Q53F19[S+80@500](NCBP3)",
        "KPNIFYSGPAS[+80]PARPR": "Q6PL18[S+80@327](ATAD2)",
        "TEFLDLDNSPLSPPS[+80]PR": "Q8NCF5[S+80@204](NFATC2IP)",
        "QGSGRES[+80]PSLASR": "Q8WWM7[S+80@339](ATXN2L)",
        # "TQLWASEPGT[+80]PPLPTSLPSQNPILK": "Q8WWM7[S+80@339](SRRM2)",
        "LQS[+80]EPESIR": "P09496[S+80@105](CLTA)",
        "RLIS[+80]PYKK": "O14929[S+80@361](HAT1)",
        "LLEDS[+80]EESSEETVSR": "O60231[S+80@103](DHX16)",
        "S[+80]PPAPGLQPM[+16]R": "P15408[S+80@200][M+16@209](FOSL2)",
        "RRLS[+80]SLR": "P62753[S+80@235](RPS6)",
        "RLS[+80]ESQLSFRR": "Q96PK6[S+80@618](RBM14)",
        "RLS[+80]LPGLLSQVSPR": "Q96Q42[S+80@483](ALS2)",
        "SPDKPGGS[+80]PSASRR": "Q9Y3T9[S+80@56](NOC2L)",
        "HLPS[+80]PPTLDSIITEYLR": "Q9Y4B6[S+80@1000](DCAF1)",
        "ST[+80]FHAGQLR": "Q7KZI7[T+80@596](MARK2)",
        "S[+80]LTNSHLEKK": "Q9H2H9[S+80@52](SLC38A1)",
        "LQTPNT[+80]FPKR": "Q14978[T+80@610](NOLC1)",
        "QIT[+80]M[+16]EELVR": "Q15149[M+16@4031][T+80@4030](PLEC)",

        "T[+56]K[+56]QTAR": "P68431[me0K@5](HIST1H3A)",
        "T[+56]K[+70]QTAR": "P68431[meK@5](HIST1H3A)",
        "T[+56]K[+28]QTAR": "P68431[me2K@5](HIST1H3A)",
        "T[+56]K[+42]QTAR-me3K": "P68431[me3K@5](HIST1H3A)",
        "T[+56]K[+42]QTAR-aK": "P68431[aK@5](HIST1H3A)",
        "K[+112.1]STGGK[+56]APR": "P68431[me0K@10][a0K@15](HIST1H3A)",
        "K[+126.1]STGGK[+56]APR": "P68431[meK@10](HIST1H3A)",
        "K[+84.1]STGGK[+56]APR": "P68431[me2K@10](HIST1H3A)",
        "K[+98.1]STGGK[+56]APR": "P68431[me3K@10](HIST1H3A)",
        "K[+98]STGGK[+56]APR": "P68431[aK@10](HIST1H3A)",
        "K[+112.1]STGGK[+42]APR": "P68431[aK@15](HIST1H3A)",
        "K[+126.1]STGGK[+42]APR": "P68431[meK@10][aK@15](HIST1H3A)",
        "K[+84.1]STGGK[+42]APR": "P68431[me2K@10][aK@15](HIST1H3A)",
        "K[+98.1]STGGK[+42]APR": "P68431[me3K@10][aK@15](HIST1H3A)",
        "K[+98]STGGK[+42]APR": "P68431[aK@10][aK@15](HIST1H3A)",
        "K[+112.1]S[+80]TGGK[+56]APR": "P68431[pS@11](HIST1H3A)",
        "K[+126.1]S[+80]TGGK[+56]APR": "P68431[meK@10][pS@11](HIST1H3A)",
        "K[+84.1]S[+80]TGGK[+56]APR": "P68431[me2K@10][pS@11](HIST1H3A)",
        "K[+98.1]S[+80]TGGK[+56]APR": "P68431[me3K@10][pS@11](HIST1H3A)",
        "K[+98]S[+80]TGGK[+56]APR": "P68431[aK@10][pS@11](HIST1H3A)",
        "K[+112.1]S[+80]TGGK[+42]APR": "P68431[pS@11][aK@15](HIST1H3A)",
        "K[+126.1]S[+80]TGGK[+42]APR": "P68431[meK@10][pS@11][aK@15](HIST1H3A)",
        "K[+84.1]S[+80]TGGK[+42]APR": "P68431[me2K@10][pS@11][aK@15](HIST1H3A)",
        "K[+98.1]S[+80]TGGK[+42]APR": "P68431[me3K@10][pS@11][aK@15](HIST1H3A)",
        "K[+98]S[+80]TGGK[+42]APR": "P68431[aK@10][pS@11][aK@15](HIST1H3A)",

        "K[+112.1]QLATK[+56]AAR": "P68431[a0K@19][a0K@24](HIST1H3A)",
        "K[+98]QLATK[+56]AAR": "P68431[aK@19](HIST1H3A)",
        "K[+112.1]QLATK[+42]AAR": "P68431[aK@24](HIST1H3A)",
        "K[+98]QLATK[+42]AAR": "P68431[aK@19][aK@24](HIST1H3A)",

        "K[+226.1]QLATK[+56]AAR": "P68431[ubK@19][a0K@24](HIST1H3A)",
        "K[+112.1]QLATK[+170.1]AAR": "P68431[ubK@24](HIST1H3A)",
        "K[+112.1]SAPATGGVK[+56]K[+56]PHR": "P68431[me0K@28][me0K@37](HIST1H3A)",
        "K[+126.1]SAPATGGVK[+56]K[+56]PHR": "P68431[meK@28][me0K@37](HIST1H3A)",
        "K[+126.1]SAPATGGVK[+70]K[+56]PHR": "P68431[meK@28][meK@37](HIST1H3A)",
        "K[+126.1]SAPATGGVK[+28]K[+56]PHR": "P68431[meK@28][me2K@37](HIST1H3A)",
        "K[+126.1]SAPATGGVK[+42]K[+56]PHR": "P68431[meK@28][me3K@37](HIST1H3A)",
        "K[+84.1]SAPATGGVK[+56]K[+56]PHR": "P68431[me2K@28][me0K@37](HIST1H3A)",
        "K[+84.1]SAPATGGVK[+70]K[+56]PHR": "P68431[me2K@28][meK@37](HIST1H3A)",
        "K[+84.1]SAPATGGVK[+28]K[+56]PHR": "P68431[me2K@28][me2K@37](HIST1H3A)",
        "K[+84.1]SAPATGGVK[+42]K[+56]PHR": "P68431[me2K@28][me3K@37](HIST1H3A)",
        "K[+98.1]SAPATGGVK[+56]K[+56]PHR": "P68431[me3K@28][me0K@37](HIST1H3A)",
        "K[+98.1]SAPATGGVK[+70]K[+56]PHR": "P68431[me3K@28][meK@37](HIST1H3A)",
        "K[+98.1]SAPATGGVK[+28]K[+56]PHR": "P68431[me3K@28][me2K@37](HIST1H3A)",
        "K[+98.1]SAPATGGVK[+42]K[+56]PHR": "P68431[me3K@28][me3K@37](HIST1H3A)",
        "K[+98]SAPATGGVK[+56]K[+56]PHR": "P68431[aK@28][me0K@37](HIST1H3A)",
        "K[+98]SAPATGGVK[+70]K[+56]PHR": "P68431[aK@28][meK@37](HIST1H3A)",
        "K[+98]SAPATGGVK[+28]K[+56]PHR": "P68431[aK@28][me2K@37](HIST1H3A)",
        "K[+98]SAPATGGVK[+42]K[+56]PHR": "P68431[aK@28][me3K@37](HIST1H3A)",
        "K[+112.1]SAPSTGGVK[+56]K[+56]PHR": "P84243[me0K@28][me0K@37](H3F3A)",
        "Y[+56]RPGTVALR": "P68431-NORM(HIST1H3A)",
        "Y[+56]QK[+56]STELLIR": "P68431[me0K@57](HIST1H3A)",
        "E[+56]IAQDFK[+56]TDLR": "P68431[me0K@80](HIST1H3A)",
        "E[+56]IAQDFK[+70]TDLR": "P68431[meK@80](HIST1H3A)",
        "E[+56]IAQDFK[+28]TDLR": "P68431[me2K@80](HIST1H3A)",
        "K[+112.1]SAPATGGVK[+70]K[+56]PHR": "P68431[me0K@28][meK@37](HIST1H3A)",
        "K[+112.1]SAPATGGVK[+28]K[+56]PHR": "P68431[me0K@28][me2K@37](HIST1H3A)",
        "K[+112.1]SAPATGGVK[+42]K[+56]PHR": "P68431[me0K@28][me3K@37](HIST1H3A)",
        "Y[+56]QK[+42]STELLIR": "P68431[aK@57](HIST1H3A)",
        "G[+56]K[+56]GGK[+56]GLGK[+56]GGAK[+56]R": "P62805[a0K@6][a0K@13][a0K@17](HIST1H4A)",
        "G[+56]K[+42]GGK[+56]GLGK[+56]GGAK[+56]R": "P62805[aK@6](HIST1H4A)",
        "G[+56]K[+56]GGK[+56]GLGK[+42]GGAK[+56]R": "P62805[aK@13](HIST1H4A)",
        "G[+56]K[+56]GGK[+56]GLGK[+56]GGAK[+42]R": "P62805[aK@17](HIST1H4A)",
        "G[+56]K[+56]GGK[+42]GLGK[+42]GGAK[+56]R": "P62805[aK@9][aK@13](HIST1H4A)",
        "G[+56]K[+42]GGK[+42]GLGK[+56]GGAK[+56]R": "P62805[aK@6][aK@9](HIST1H4A)",
        "G[+56]K[+42]GGK[+56]GLGK[+56]GGAK[+42]R": "P62805[aK@6][aK@17](HIST1H4A)",
        "G[+56]K[+56]GGK[+56]GLGK[+42]GGAK[+42]R": "P62805[aK@13][aK@17](HIST1H4A)",
        "G[+56]K[+56]GGK[+42]GLGK[+42]GGAK[+42]R": "P62805[aK@9][aK@13][aK@17](HIST1H4A)",
        "G[+56]K[+42]GGK[+42]GLGK[+42]GGAK[+56]R": "P62805[aK@6][aK@9][aK@13](HIST1H4A)",
        "G[+56]K[+42]GGK[+42]GLGK[+56]GGAK[+42]R": "P62805[aK@6][aK@9][aK@17](HIST1H4A)",
        "G[+56]K[+42]GGK[+42]GLGK[+42]GGAK[+42]R": "P62805[aK@6][aK@9][aK@13][aK@17](HIST1H4A)",
        "G[+56]K[+56]GGK[+56]GLGK[+56]GGAK[+70]R": "P62805[meK@17](HIST1H4A)",
        "K[+112.1]VLR": "P62805[me0K@21](HIST1H4A)",
        "K[+126.1]VLR": "P62805[meK@21](HIST1H4A)",
        "K[+84.1]VLR": "P62805[me2K@21](HIST1H4A)",
        "K[+98.1]VLR": "P62805[me3K@21](HIST1H4A)",
        "D[+56]AVTYTEHAK[+56]R": "P62805-NORM(HIST1H4A)",
        "Y[+56]QK[+28]STELLIR": "P68431[me2K@57](HIST1H3A)"
    }
    P100_data = pd.read_csv('static/data/P100-all-plates-Level4.2018_08_07.processed.csv', sep=',', header=None)
    p100_unique_complete = {}
    for column in range(21, P100_data.shape[1]):
        # print (column)
        # identifier = str(P100_data[column][10])  +"++"+str(P100_data[column][3])  +"++"+ str(P100_data[column][7]) +"++"+ str(P100_data[column][11])
        p100_identifier_complete = str(P100_data[column][1]) + "++" + str(P100_data[column][3]) + "++" + str(
            P100_data[column][4]) + "++"  + str(P100_data[column][7]) + "++" + str(P100_data[column][8]) + "++" + \
                                   str(P100_data[column][11]) + "++" + str(P100_data[column][12]) + "++" + \
                                   str(P100_data[column][13]) + "++" + str(P100_data[column][14]) + "++" \
                                   + str(P100_data[column][15]) + "++" + \
                                   str(P100_data[column][16]) + "++" + str(P100_data[column][17]) + "++" + \
                                   str(P100_data[column][18]) + "++" + str(P100_data[column][19]) + "++" + \
                                   str(P100_data[column][20]) + "++" + str(P100_data[column][21]) + "++" + \
                                   str(P100_data[column][22]) + "++" + str(P100_data[column][23]) + "++" + \
                                    str(P100_data[column][24]) + "++" + str(P100_data[column][25]) + "++" + \
                                    str(P100_data[column][26]) + "++" + str(P100_data[column][27]) + "++" + \
                                    str(P100_data[column][28])

        if p100_identifier_complete not in p100_unique_complete.keys():
            # print identifier_complete
            p100_unique_complete[p100_identifier_complete] = []
            p100_unique_complete.get(p100_identifier_complete).append(column)
        else:
            p100_unique_complete.get(p100_identifier_complete).append(column)


    iterator = 0
    for key in p100_unique_complete:
        val = p100_unique_complete.get(key)
        value = list(val)
        # print key
        for item in value:
            iterator += 1
        #     print item
        # print "------------------------"

    print iterator
    #print unique_complete_modified

    total = 0
    total_complete = 0


    print "size of unique_complete keys"
    print len(p100_unique_complete.keys())
    print "size of unique_complete"
    print len(p100_unique_complete)

    unique_sum = {}
    # unique_sum_modified = {}
    unique_pvalue = {}


    col_cell_line = []
    col_pert_name = []
    col_pert_dose = []
    col_pert_time = []
    col_pert_type = []
    col_pert_signature = []

    my_data_con = []
    keyIter = 0
    first_loop = True

    #{"pertAnnotation": [{"LINCS_UNIQUE": "LINCSTP_0094", "pertVehicle": "DMSO", "pertName": "CPI-169", "cellId": "A375",
     #                    "pubchemCid": "71712226", "pertId": "BRD-A48881734", "pertType": "trt_cp", "pertTime": "3",
      #                   "pertDose": "1", "id": "P1", "lsmId": "LSM-42750"},
       #                 {"LINCS_UNIQUE": "LINCSTP_0236", "pertVehicle": "DMSO", "pertName": "CPI-169", "cellId": "MCF7",
        #                 "pubchemCid": "71712226", "pertId": "BRD-A48881734", "pertType": "trt_cp", "pertTime": "3",
         #                "pertDose": "1", "id": "P2", "lsmId": "LSM-42750"},
          #              {"LINCS_UNIQUE": "LINCSTP_0259", "pertVehicle": "DMSO", "pertName": "CPI-169", "cellId": "NPC",
           #              "pubchemCid": "71712226", "pertId": "BRD-A48881734", "pertType": "trt_cp", "pertTime": "3",
            #             "pertDose": "1", "id": "P3", "lsmId": "LSM-42750"},
             #           {"LINCS_UNIQUE": "LINCSTP_0116", "pertVehicle": "DMSO", "pertName": "CPI-169", "cellId": "A549",
              #           "pubchemCid": "71712226", "pertId": "BRD-A48881734", "pertType": "trt_cp", "pertTime": "3",
               #          "pertDose": "1", "id": "P4", "lsmId": "LSM-42750"},
                #        {"LINCS_UNIQUE": "LINCSTP_0275", "pertVehicle": "DMSO", "pertName": "CPI-169", "cellId": "PC3",
                 #        "pubchemCid": "71712226", "pertId": "BRD-A48881734", "pertType": "trt_cp", "pertTime": "3",
                  #       "pertDose": "1", "id": "P5", "lsmId": "LSM-42750"},
                   #     {"LINCS_UNIQUE": "LINCSTP_0441", "pertVehicle": "DMSO", "pertName": "CPI-169", "cellId": "YAPC",
                    #     "pubchemCid": "71712226", "pertId": "BRD-A48881734", "pertType": "trt_cp", "pertTime": "3",
                         # "pertDose": "1", "id": "P6", "lsmId": "LSM-42750"}]}, "LINCSTP_0002":

    pertAnnotation = {}
    # this is for slicing the clustergrammer
    pertAnnotation_for_slicing = {}
    pertJson = {}
    unique_iterator = 1
    keylist = p100_unique_complete.keys()
    keylist.sort()
    for key in keylist:
        initial = [0.0]*96
        initial2 = [0.0] * 96
        initial_not_na = [0.0] * 96
        signature_id = "LINCS_TP" + str(unique_iterator)
        unique_iterator += 1
        unique_sum[key] = initial
        unique_pvalue[key] = initial2
        val = p100_unique_complete.get(key)
        key_splitted = key.split("++")
        print key_splitted
        col_cell_line.append('Cell Line: ' + key_splitted[0])
        col_pert_dose.append('Dose: ' + key_splitted[5] + key_splitted[6])
        col_pert_time.append('Time: ' + key_splitted[9] + key_splitted[10])
        col_pert_name.append('Perturbations: ' + key_splitted[8]+ "/"+key_splitted[5] + key_splitted[6]+ "/"+key_splitted[9] + key_splitted[10])
        col_pert_signature.append('Signature ID: ' + signature_id)
        if "_sg01" not in key_splitted[8] and "_sg02" not in key_splitted[8]:
            cpType = 'Small Molecule'
        else:
            cpType = 'CRISPR'

        col_pert_type.append('Type: ' + cpType)
        compound_name = str(key_splitted[8]).upper()
        cp_metaData = {}
        cp_metaData_for_data_slicing = {}
        cp_metaData["Cell_line"] = key_splitted[0]
        cp_metaData["Dose"] = key_splitted[5] + " "+ key_splitted[6]
        cp_metaData["Time"] = key_splitted[9] + " "+ key_splitted[10]
        cp_metaData["Perturbations"] = compound_name
        cp_metaData["Type"] = cpType
        cp_metaData["Sig_id"] = signature_id

        cp_metaData_for_data_slicing["Cell_line"] = key_splitted[0]
        cp_metaData_for_data_slicing["Dose"] = key_splitted[5] + " "+ key_splitted[6]
        cp_metaData_for_data_slicing["Time"] = key_splitted[9] + " "+ key_splitted[10]
        cp_metaData_for_data_slicing["Perturbations"] = compound_name
        cp_metaData_for_data_slicing["Type"] = cpType
        cp_metaData_for_data_slicing["Sig_id"] = signature_id

        if compound_name not in pertAnnotation.keys():
            # print identifier_complete
            pertAnnotation[compound_name] = []
            pertAnnotation.get(compound_name).append(cp_metaData)
        else:
            pertAnnotation.get(compound_name).append(cp_metaData)

        item = list(val)
        # if len(item) == 1:
        #     print '111111   %d  =========== %s %s' % (keyIter, key, val)

        keyIter += 1

        value_list = 0
        sig2PeptideValue = {}
        for item_iter in range(0, 96):

            aux_list = []
            for col in item:

                if str(P100_data[col][item_iter + 29]) != "nan":
                    aux_list.append(float(P100_data[col][item_iter + 29]))
                    unique_sum.get(key)[item_iter] += float(P100_data[col][item_iter + 29])
                    initial_not_na[item_iter] += 1
                else:
                    print cp_metaData
            if initial_not_na[item_iter] > 0:
                unique_sum.get(key)[item_iter] /= initial_not_na[item_iter]




            shorthand = petide2ptmAndGene[str(P100_data[6][29 + item_iter])]
            sig2PeptideValue[shorthand] = unique_sum.get(key)[item_iter]

            meanValue = unique_sum.get(key)[item_iter]
            average = unique_sum.get(key)[item_iter]

            if len(aux_list) > 1:
                variance2 = sum((average - value) ** 2 for value in aux_list) / (len(aux_list) - 1)
                t_stat2 = meanValue / np.sqrt(variance2 / float(len(aux_list)))
                # p_val = stats.t.sf(abs(t_stat), len(item)-1)*2
                p_val2 = stats.t.sf(abs(t_stat2), len(aux_list) - 1) * 2
            else:
                p_val2 = None
            unique_pvalue.get(key)[item_iter] = p_val2
            # variance = np.var(aux_list)

            # print aux_list
            # t_stat = meanValue/np.sqrt(variance/float(len(item)))













            if (average != unique_sum.get(key)[item_iter] or average != meanValue or meanValue != unique_sum.get(key)[item_iter]):
                print average
                print meanValue
                print unique_sum.get(key)[item_iter]
                print "++++++++++++++++++"
            if len(aux_list) > 1:

                variance2 = sum((average - value) ** 2 for value in aux_list) / (len(aux_list)-1)
        #print "end of loop"
        cp_metaData_for_data_slicing["data"] = unique_sum.get(key)

        if compound_name not in pertAnnotation_for_slicing.keys():
            # print identifier_complete
            pertAnnotation_for_slicing[compound_name] = []
            pertAnnotation_for_slicing.get(compound_name).append(cp_metaData_for_data_slicing)
        else:
            pertAnnotation_for_slicing.get(compound_name).append(cp_metaData_for_data_slicing)


        x1 = np.array(unique_sum.get(key))
        #print "end of loop1"
        x2 = x1[:, np.newaxis]
        #print "end of loop2"
        my_data_final = np.asarray(x2)
        #print "end of loop3"
        if first_loop:
            my_data_con = my_data_final
            first_loop = False

        else:
            my_data_con = np.hstack((my_data_con, my_data_final))
        print "end of loop over unique"
        pertJson[signature_id] = sig2PeptideValue

    for keys in pertAnnotation:
        pertJson[keys] = pertAnnotation[keys]





    print "======================"
    print "======================"
    print "======================"
    print my_data_con.shape




    row_peptide = []
    row_protein = []
    row_gene = []
    row_shorthand = []

    for i in range(0, 96):
        row_peptide.append('Peptide: ' + str(P100_data[6][29 + i]))
        shorthand = petide2ptmAndGene[str(P100_data[6][29 + i])]
        row_shorthand.append('PTM Proteins: ' + shorthand)
        row_protein.append('Protein: ' + str(P100_data[11][29 + i]))
        row_gene.append('Gene: ' + str(P100_data[2][29 + i]))


    arrays_rows2 = [np.array(row_shorthand),np.array(row_peptide),np.array(row_protein),np.array(row_gene)]
    arrays_columns2 = [np.array(col_pert_name), np.array(col_pert_signature), np.array(col_pert_time), np.array(col_pert_dose), np.array(col_pert_type),np.array(col_cell_line)]

    np.save("static/data/P100_rowMetaDat", arrays_rows2)

    tuples_rows2 = list(zip(*arrays_rows2))
    #print tuples_rows2
    tuples_columns2 = list(zip(*arrays_columns2))
    #print tuples_columns2
    rows_labels2 = pd.MultiIndex.from_tuples(tuples_rows2)
    columns_labels2 = pd.MultiIndex.from_tuples(tuples_columns2)

    #my_data_con = my_data_con.fillna(0)
    print "after na to zero"
    a2 = np.nan_to_num(my_data_con)
    #a2 = my_data_con
    aa = a2.astype(np.float)
    print a2
    print a2.shape
    print columns_labels2.shape
    print rows_labels2.shape
    # {'nop': row1, 'o0p': row2, 'zaz': row3, 'zax': row4, 'oof': row5, 'oye': row6}
    #df = pd.DataFrame(a, index=rows_labels, columns=columns_labels)
    df2 = pd.DataFrame(aa, index=rows_labels2, columns=columns_labels2)
    #print df2
    P100_all_json = loadFile.make_json_from_txt(df2)
    json.dump(P100_all_json, open("static/data/P100_aggregated_clustergram.json", 'w'))
    json.dump(pertJson, open("static/data/P100_processed_perturb_new.json", 'w'))
    json.dump(pertAnnotation_for_slicing, open("static/data/P100_processed_perturb_for_clustergramm_slicing.json", 'w'))

    # str(P100_data[column][1]) + "++" + str(P100_data[column][3]) + "++" + str(
    #     P100_data[column][4]) + "++" + str(P100_data[column][7]) + "++" + \
    # str(P100_data[column][8]) + "++" + \
    # str(P100_data[column][11]) + "++" + str(P100_data[column][12]) + "++" + \
    # str(P100_data[column][13]) + "++" + str(P100_data[column][14]) + "++" \
    # + str(P100_data[column][15]) + "++" + \
    # str(P100_data[column][16]) + "++" + str(P100_data[column][17]) + "++" + \
    # str(P100_data[column][18]) + "++" + str(P100_data[column][19]) + "++" + \
    # str(P100_data[column][20]) + "++" + str(P100_data[column][21]) + "++" + \
    # str(P100_data[column][22]) + "++" + str(P100_data[column][23]) + "++" + \
    # str(P100_data[column][24]) + "++" + str(P100_data[column][25]) + "++" + \
    # str(P100_data[column][26]) + "++" + str(P100_data[column][27]) + "++" + \
    # str(P100_data[column][28])

    #making p-value file
    first_col = []
    first_col.append("")
    first_col.append("LINCS_Signature_ID")
    first_col.append(str(P100_data[0][1]))
    first_col.append(str(P100_data[0][3]))
    first_col.append(str(P100_data[0][4]))
    first_col.append(str(P100_data[0][7]))
    first_col.append(str(P100_data[0][8]))
    first_col.append(str(P100_data[0][11]))
    first_col.append(str(P100_data[0][12]))
    first_col.append(str(P100_data[0][13]))
    first_col.append(str(P100_data[0][14]))
    first_col.append(str(P100_data[0][15]))
    first_col.append(str(P100_data[0][16]))
    first_col.append(str(P100_data[0][17]))
    first_col.append(str(P100_data[0][18]))
    first_col.append(str(P100_data[0][19]))
    first_col.append(str(P100_data[0][20]))
    first_col.append(str(P100_data[0][21]))
    first_col.append(str(P100_data[0][22]))
    first_col.append(str(P100_data[0][23]))
    first_col.append(str(P100_data[0][24]))
    first_col.append(str(P100_data[0][25]))
    first_col.append(str(P100_data[0][26]))
    first_col.append(str(P100_data[0][27]))
    first_col.append(str(P100_data[0][28]))

    for i in range(0, 96):
        first_col.append(str(P100_data[0][29 + i]))

    x1 = np.array(first_col)
    x2 = x1[:, np.newaxis]
    my_data_final = np.asarray(x2)

    for col in range(1, 11):
        ith_col = []
        ith_col.append(str(P100_data[col][0]))
        ith_col.append(str(""))
        ith_col.append(str(""))
        ith_col.append(str(""))
        ith_col.append(str(""))
        ith_col.append(str(""))
        ith_col.append(str(""))
        ith_col.append(str(""))
        ith_col.append(str(""))
        ith_col.append(str(""))
        ith_col.append(str(""))
        ith_col.append(str(""))
        ith_col.append(str(""))
        ith_col.append(str(""))
        ith_col.append(str(""))
        ith_col.append(str(""))
        ith_col.append(str(""))
        ith_col.append(str(""))
        ith_col.append(str(""))
        ith_col.append(str(""))
        ith_col.append(str(""))
        ith_col.append(str(""))
        ith_col.append(str(""))
        ith_col.append(str(""))
        ith_col.append(str(""))

        for i in range(0, 96):
            ith_col.append(str(P100_data[col][29 + i]))

        x1 = np.array(ith_col)
        x2 = x1[:, np.newaxis]
        col_data_final = np.asarray(x2)
        my_data_final = np.hstack((my_data_final, col_data_final))

    #Adding PTMProtein
    ith_col = []
    ith_col.append(str("ShortHand"))
    ith_col.append(str(""))
    ith_col.append(str(""))
    ith_col.append(str(""))
    ith_col.append(str(""))
    ith_col.append(str(""))
    ith_col.append(str(""))
    ith_col.append(str(""))
    ith_col.append(str(""))
    ith_col.append(str(""))
    ith_col.append(str(""))
    ith_col.append(str(""))
    ith_col.append(str(""))
    ith_col.append(str(""))
    ith_col.append(str(""))
    ith_col.append(str(""))
    ith_col.append(str(""))
    ith_col.append(str(""))
    ith_col.append(str(""))
    ith_col.append(str(""))
    ith_col.append(str(""))
    ith_col.append(str(""))
    ith_col.append(str(""))
    ith_col.append(str(""))
    ith_col.append(str(""))

    for i in range(0, 96):

        shorthand = petide2ptmAndGene[str(P100_data[6][29 + i])]
        ith_col.append(shorthand)
    x1 = np.array(ith_col)
    x2 = x1[:, np.newaxis]
    col_data_final = np.asarray(x2)
    my_data_final = np.hstack((my_data_final, col_data_final))
    my_data_final2 = my_data_final
    my_data_final3 = my_data_final
    my_pvalue_con = my_data_final2
    my_average_con = my_data_final3
    print "my_data_final.shape"
    print my_data_final.shape
    unique_iterator = 1
    keylist = p100_unique_complete.keys()
    keylist.sort()
    for key in keylist:


        signature_id = "LINCS_TP" + str(unique_iterator)
        unique_iterator += 1

        key_splitted = key.split("++")

        key_list1 = []
        key_list1.append("")
        key_list1.append(signature_id)



        for list_member in key_splitted:
            key_list1.append(list_member)

        # change list to array
        key_arr1 = np.array(key_list1)

        # This is for making data n by 1 instead of 1 by n
        key_arr2 = key_arr1[:, np.newaxis]


        x1 = np.array(unique_sum.get(key))
        x2 = x1[:, np.newaxis]

        my_data = np.asarray(x2)
        col_data = np.vstack((key_arr2, my_data))

        # Repeating for p-values
        x3 = np.array(unique_pvalue.get(key))
        x4 = x3[:, np.newaxis]

        my_data_pvalue = np.asarray(x4)
        col_data_pvalue = np.vstack((key_arr2, my_data_pvalue))

        my_average_con = np.hstack((my_average_con, col_data))
        my_pvalue_con = np.hstack((my_pvalue_con, col_data_pvalue))














    print "my_data.shape"
    print my_pvalue_con.shape

    np.savetxt("static/data/P100_average.txt", my_average_con,delimiter='\t', fmt="%s")
    np.savetxt("static/data/P100_pValue.txt", my_pvalue_con,delimiter='\t', fmt="%s")







    return P100_all_json

@app.route("/api/printptm")
def print_ptm():
    peptide2ptm = {
    # "IYQY[+80]IQSR": "Q13627[Y+80@321]",
    # "TPKDS[+80]PGIPPSANAHQLFR": "P51812[S+80@369]",
    # "RNS[+80]SEASSGDFLDLK": "Q9UK76[S+80@87]",
    # "LPLVPES[+80]PRR": "Q86WB0[S+80@321]",
    # "ANAS[+80]PQKPLDLK": "Q9Y618[S+80@956]",
    # "LENS[+80]PLGEALR": "Q9NX40[S+80@108]",
    # "ANS[+80]FVGTAQYVSPELLTEK": "O15530[S+80@241]",
    # "TNPPTQKPPS[+80]PPMSGR": "Q8IZP0[S+80@183]",
    # "SNS[+80]LPHSAVSNAGSK": "Q8TBZ3[S+80@434]",
    # "VGS[+80]LDNVGHLPAGGAVK": "P27816[S+80@1073]",
    # "AAPEAS[+80]SPPASPLQHLLPGK": "Q96TA1[S+80@691]",
    # "S[+122]DKPDM[+16]AEIEKFDK": "P62328[S+122@2][M+16@7]",
    # "S[+122]DKPDMAEIEKFDK": "P62328[S+122@2]",
    # "SLS[+80]LGDKEISR": "Q9UMZ2[S+80@1075]",
    # "DLVQPDKPAS[+80]PK": "Q6PJT7[S+80@515]",
    # "SPS[+80]PAHLPDDPKVAEK": "Q92615[S+80@601]",
    # "S[+80]IQDLTVTGTEPGQVSSR": "O43318[S+80@439]",
    # "IHS[+80]PIIR": "O60885[S+80@1117]",
    # "TFS[+80]LTEVR": "O95239[S+80@801]",
    # "SLVGS[+80]WLK": "Q6ICG6[S+80@362]",
    # "S[+80]PPAPGLQPMR": "P15408[S+80@200]",
    # "LAS[+80]PELER": "P17535[S+80@100]",
    # "IGPLGLS[+80]PK": "P30050[S+80@38]",
    # "TPS[+80]IQPSLLPHAAPFAK": "P35658[S+80@1023]",
    # "HAS[+80]PILPITEFSDIPR": "P42167[S+80@306]",
    # "LIPGPLS[+80]PVAR": "P48634[S+80@1219]",
    # "LGM[+16]LS[+80]PEGTC[+57]K": "P49327[S+80@207][M+16@205][C+57@212]",
    # "LGMLS[+80]PEGTC[+57]K": "P49327[S+80@207][C+57@212]",
    # "ISNLS[+80]PEEEQGLWK": "Q5HYJ3[S+80@193]",
    # "VSMPDVELNLKS[+80]PK": "Q09666[S+80@3426]",
    # "S[+122]DNGELEDKPPAPPVR": "Q13177[S+122@2]",
    # "KAYS[+80]FC[+57]GTVEYM[+16]APEVVNR": "Q15418[S+80@221][M+16@229][C+57@223]",
    # "KAYS[+80]FC[+57]GTVEYMAPEVVNR": "Q15418[S+80@221][C+57@223]",
    # "NDS[+80]WGSFDLR": "Q7Z417[S+80@652]",
    # "LEVTEIVKPS[+80]PK": "Q7Z6E9[S+80@1179]",
    # "YGS[+80]PPQRDPNWNGER": "O15234[S+80@265]",
    # "QDDS[+80]PPRPIIGPALPPGFIK": "Q8IXQ4[S+80@105]",
    # "SFS[+80]ADNFIGIQR": "Q8N7R7[S+80@344]",
    # "VLS[+80]PLIIK": "Q8NCN4[S+80@403]",
    # "AGS[+80]PDVLR": "Q8NDX6[S+80@44]",
    # "LGPGRPLPTFPTSEC[+57]TS[+80]DVEPDTR": "Q8TDD1[S+80@75][C+57@73]",
    # "LAAPSVSHVS[+80]PR": "Q8WXE1[S+80@224]",
    # "VDDDS[+80]LGEFPVTNSR": "Q92785[S+80@142]",
    # "NEEPVRS[+80]PERR": "Q92922[S+80@310]",
    # "LFIIRGS[+80]PQQIDHAK": "Q92945[S+80@480]",
    # "S[+80]IEVENDFLPVEK": "Q96B97[S+80@230]",
    # "TAPTLS[+80]PEHWK": "Q96JM3[S+80@405]",
    # "VLS[+80]PTAAKPSPFEGK": "Q96QC0[S+80@313]",
    # "SSDQPLTVPVS[+80]PK": "Q9ULW0[S+80@738]",
    # "FYETKEESYS[+80]PSKDR": "Q96T23[S+80@473]",
    # "SDS[+80]PENKYSDSTGHSK": "Q9BTA9[S+80@64]",
    # "S[+80]IPLSIK": "Q9C0C9[S+80@515]",
    # "RLS[+80]QSDEDVIR": "Q9H7D7[S+80@121]",
    # "ATS[+80]PVKSTTSITDAK": "Q9NQW6[S+80@295]",
    # "ALGS[+80]PTKQLLPC[+57]EMAC[+57]NEK": "Q9NR45[S+80@275][C+57@283][C+57@287]",
    # "YLLGDAPVS[+80]PSSQK": "Q9NYB0[S+80@203]",
    # "ANS[+80]PEKPPEAGAAHKPR": "Q9UFC0[S+80@212]",
    # "SEVQQPVHPKPLS[+80]PDSR": "Q9UHB6[S+80@362]",
    # "ETPHS[+80]PGVEDAPIAK": "Q9UHB6[S+80@490]",
    # "SQS[+80]PHYFR": "Q9UKJ3[S+80@1035]",
    # "DRS[+80]SPPPGYIPDELHQVAR": "Q9Y2U5[S+80@163]",
    # "SPALKS[+80]PLQSVVVR": "Q9Y2W1[S+80@253]",
    # "AFGSGIDIKPGT[+80]PPIAGR": "Q9Y520[T+80@2673]",
    # "SFS[+80]SQRPVDR": "Q9Y520[S+80@1544]",
    # "VYT[+80]HEVVTLWYR": "P06493[T+80@161]",
    # "SST[+80]PLPTISSSAENTR": "P42167[T+80@160]",
    # "QIT[+80]MEELVR": "Q15149[T+80@4030]",
    # "TQLWASEPGT[+80]PPLPTSLPSQNPILK": "Q9BXP5[T+80@544]",
    # # "TQLWASEPGT[+80]PPLPTSLPSQNPILK": "Q8WWM7[S+80@339]",
    # "ALPQT[+80]PRPR": "Q9UQ35[T+80@1492]",
    # "SMS[+80]VDLSHIPLKDPLLFK": "A0JNW5[S+80@935]",
    # "S[+80]PTGPSNSFLANMGGTVAHK": "Q96I25[S+80@222]",
    # "S[+80]LTAHSLLPLAEK": "Q86VI3[S+80@1424]",
    # "S[+80]FAGNLNTYKR": "Q01813[S+80@386]",
    # "HRPS[+80]PPATPPPK": "Q8IYB3[S+80@402]",
    # "LHS[+80]APNLSDLHVVRPK": "O75385[S+80@556]",
    # "TLGRRDS[+80]SDDWEIPDGQITVGQR": "P15056[S+80@446]",
    # "A[+42]TTATM[+16]ATSGS[+80]AR": "P38919[S+80@12][M+16@7][A+42@2]",
    # "A[+42]TTATMATSGS[+80]AR": "P38919[S+80@12][A+42@2]",
    # "IHVSRS[+80]PTRPR": "Q499Z4[S+80@189]",
    # "RPHS[+80]PEKAFSSNPVVR": "Q53F19[S+80@500]",
    # "KPNIFYSGPAS[+80]PARPR": "Q6PL18[S+80@327]",
    # "TEFLDLDNSPLSPPS[+80]PR": "Q8NCF5[S+80@204]",
    # "QGSGRES[+80]PSLASR": "Q8WWM7[S+80@339]",
    # "LQS[+80]EPESIR": "P09496[S+80@105]",
    # "RLIS[+80]PYKK": "O14929[S+80@361]",
    # "LLEDS[+80]EESSEETVSR": "O60231[S+80@103]",
    # "S[+80]PPAPGLQPM[+16]R": "P15408[S+80@200][M+16@209]",
    # "RRLS[+80]SLR": "P62753[S+80@235]",
    # "RLS[+80]ESQLSFRR": "Q96PK6[S+80@618]",
    # "RLS[+80]LPGLLSQVSPR": "Q96Q42[S+80@483]",
    # "SPDKPGGS[+80]PSASRR": "Q9Y3T9[S+80@56]",
    # "HLPS[+80]PPTLDSIITEYLR": "Q9Y4B6[S+80@1000]",
    # "ST[+80]FHAGQLR": "Q7KZI7[T+80@596]",
    # "S[+80]LTNSHLEKK": "Q9H2H9[S+80@52]",
    # "LQTPNT[+80]FPKR": "Q14978[T+80@610]",
    # "QIT[+80]M[+16]EELVR": "Q15149[M+16@4031][T+80@4030]",

    "T[+56]K[+56]QTAR": "P68431[me0K@5]",
    "T[+56]K[+70]QTAR": "P68431[meK@5]",
    "T[+56]K[+28]QTAR": "P68431[me2K@5]",
    "T[+56]K[+42]QTAR-me3K": "P68431[me3K@5]",
    "T[+56]K[+42]QTAR-aK": "P68431[aK@5]",
    "K[+112.1]STGGK[+56]APR": "P68431[me0K@10][a0K@15]",
    "K[+126.1]STGGK[+56]APR": "P68431[meK@10]",
    "K[+84.1]STGGK[+56]APR": "P68431[me2K@10]",
    "K[+98.1]STGGK[+56]APR": "P68431[me3K@10]",
    "K[+98]STGGK[+56]APR": "P68431[aK@10]",
    "K[+112.1]STGGK[+42]APR": "P68431[aK@15]",
    "K[+126.1]STGGK[+42]APR": "P68431[meK@10][aK@15]",
    "K[+84.1]STGGK[+42]APR": "P68431[me2K@10][aK@15]",
    "K[+98.1]STGGK[+42]APR": "P68431[me3K@10][aK@15]",
    "K[+98]STGGK[+42]APR": "P68431[aK@10][aK@15]",
    "K[+112.1]S[+80]TGGK[+56]APR": "P68431[pS@11]",
    "K[+126.1]S[+80]TGGK[+56]APR": "P68431[meK@10][pS@11]",
    "K[+84.1]S[+80]TGGK[+56]APR": "P68431[me2K@10][pS@11]",
    "K[+98.1]S[+80]TGGK[+56]APR": "P68431[me3K@10][pS@11]",
    "K[+98]S[+80]TGGK[+56]APR": "P68431[aK@10][pS@11]",
    "K[+112.1]S[+80]TGGK[+42]APR": "P68431[pS@11][aK@15]",
    "K[+126.1]S[+80]TGGK[+42]APR": "P68431[meK@10][pS@11][aK@15]",
    "K[+84.1]S[+80]TGGK[+42]APR": "P68431[me2K@10][pS@11][aK@15]",
    "K[+98.1]S[+80]TGGK[+42]APR": "P68431[me3K@10][pS@11][aK@15]",
    "K[+98]S[+80]TGGK[+42]APR": "P68431[aK@10][pS@11][aK@15]",

    "K[+112.1]QLATK[+56]AAR": "P68431[a0K@19][a0K@24]",
    "K[+98]QLATK[+56]AAR": "P68431[aK@19]",
    "K[+112.1]QLATK[+42]AAR": "P68431[aK@24]",
    "K[+98]QLATK[+42]AAR": "P68431[aK@19][aK@24]",

    "K[+226.1]QLATK[+56]AAR": "P68431[ubK@19][a0K@24]",
    "K[+112.1]QLATK[+170.1]AAR": "P68431[ubK@24]",
    "K[+112.1]SAPATGGVK[+56]K[+56]PHR": "P68431[me0K@28][me0K@37]",
    "K[+126.1]SAPATGGVK[+56]K[+56]PHR": "P68431[meK@28][me0K@37]",
    "K[+126.1]SAPATGGVK[+70]K[+56]PHR": "P68431[meK@28][meK@37]",
    "K[+126.1]SAPATGGVK[+28]K[+56]PHR": "P68431[meK@28][me2K@37]",
    "K[+126.1]SAPATGGVK[+42]K[+56]PHR": "P68431[meK@28][me3K@37]",
    "K[+84.1]SAPATGGVK[+56]K[+56]PHR": "P68431[me2K@28][me0K@37]",
    "K[+84.1]SAPATGGVK[+70]K[+56]PHR": "P68431[me2K@28][meK@37]",
    "K[+84.1]SAPATGGVK[+28]K[+56]PHR": "P68431[me2K@28][me2K@37]",
    "K[+84.1]SAPATGGVK[+42]K[+56]PHR": "P68431[me2K@28][me3K@37]",
    "K[+98.1]SAPATGGVK[+56]K[+56]PHR": "P68431[me3K@28][me0K@37]",
    "K[+98.1]SAPATGGVK[+70]K[+56]PHR": "P68431[me3K@28][meK@37]",
    "K[+98.1]SAPATGGVK[+28]K[+56]PHR": "P68431[me3K@28][me2K@37]",
    "K[+98.1]SAPATGGVK[+42]K[+56]PHR": "P68431[me3K@28][me3K@37]",
    "K[+98]SAPATGGVK[+56]K[+56]PHR": "P68431[aK@28][me0K@37]",
    "K[+98]SAPATGGVK[+70]K[+56]PHR": "P68431[aK@28][meK@37]",
    "K[+98]SAPATGGVK[+28]K[+56]PHR": "P68431[aK@28][me2K@37]",
    "K[+98]SAPATGGVK[+42]K[+56]PHR": "P68431[aK@28][me3K@37]",
    "K[+112.1]SAPSTGGVK[+56]K[+56]PHR": "P84243[me0K@28][me0K@37]",
    "Y[+56]RPGTVALR": "P68431-NORM",
    "Y[+56]QK[+56]STELLIR": "P68431[me0K@57]",
    "E[+56]IAQDFK[+56]TDLR": "P68431[me0K@80]",
    "E[+56]IAQDFK[+70]TDLR": "P68431[meK@80]",
    "E[+56]IAQDFK[+28]TDLR": "P68431[me2K@80]",
    "K[+112.1]SAPATGGVK[+70]K[+56]PHR": "P68431[me0K@28][meK@37]",
    "K[+112.1]SAPATGGVK[+28]K[+56]PHR": "P68431[me0K@28][me2K@37]",
    "K[+112.1]SAPATGGVK[+42]K[+56]PHR": "P68431[me0K@28][me3K@37]",
    "Y[+56]QK[+42]STELLIR": "P68431[aK@57]",
    "G[+56]K[+56]GGK[+56]GLGK[+56]GGAK[+56]R": "P62805[a0K@6][a0K@13][a0K@17]",
    "G[+56]K[+42]GGK[+56]GLGK[+56]GGAK[+56]R": "P62805[aK@6]",
    "G[+56]K[+56]GGK[+56]GLGK[+42]GGAK[+56]R": "P62805[aK@13]",
    "G[+56]K[+56]GGK[+56]GLGK[+56]GGAK[+42]R": "P62805[aK@17]",
    "G[+56]K[+56]GGK[+42]GLGK[+42]GGAK[+56]R": "P62805[aK@9][aK@13]",
    "G[+56]K[+42]GGK[+42]GLGK[+56]GGAK[+56]R": "P62805[aK@6][aK@9]",
    "G[+56]K[+42]GGK[+56]GLGK[+56]GGAK[+42]R": "P62805[aK@6][aK@17]",
    "G[+56]K[+56]GGK[+56]GLGK[+42]GGAK[+42]R": "P62805[aK@13][aK@17]",
    "G[+56]K[+56]GGK[+42]GLGK[+42]GGAK[+42]R": "P62805[aK@9][aK@13][aK@17]",
    "G[+56]K[+42]GGK[+42]GLGK[+42]GGAK[+56]R": "P62805[aK@6][aK@9][aK@13]",
    "G[+56]K[+42]GGK[+42]GLGK[+56]GGAK[+42]R": "P62805[aK@6][aK@9][aK@17]",
    "G[+56]K[+42]GGK[+42]GLGK[+42]GGAK[+42]R": "P62805[aK@6][aK@9][aK@13][aK@17]",
    "G[+56]K[+56]GGK[+56]GLGK[+56]GGAK[+70]R": "P62805[meK@17]",
    "K[+112.1]VLR": "P62805[me0K@21]",
    "K[+126.1]VLR": "P62805[meK@21]",
    "K[+84.1]VLR": "P62805[me2K@21]",
    "K[+98.1]VLR": "P62805[me3K@21]",
    "D[+56]AVTYTEHAK[+56]R": "P62805-NORM",
    "Y[+56]QK[+28]STELLIR": "P68431[me2K@57]"

    }
    petide2ptmAndGene = {
        # "IYQY[+80]IQSR": "Q13627[Y+80@321](DYRK1A)",
        # "TPKDS[+80]PGIPPSANAHQLFR": "P51812[S+80@369](RPS6KA3)",
        # "RNS[+80]SEASSGDFLDLK": "Q9UK76[S+80@87](JPT1)",
        # "LPLVPES[+80]PRR": "Q86WB0[S+80@321](ZC3HC1)",
        # "ANAS[+80]PQKPLDLK": "Q9Y618[S+80@956](NCOR2)",
        # "LENS[+80]PLGEALR": "Q9NX40[S+80@108](OCIAD1)",
        # "ANS[+80]FVGTAQYVSPELLTEK": "O15530[S+80@241](PDPK1)",
        # "TNPPTQKPPS[+80]PPMSGR": "Q8IZP0[S+80@183](ABI1)",
        # "SNS[+80]LPHSAVSNAGSK": "Q8TBZ3[S+80@434](WDR20)",
        # "VGS[+80]LDNVGHLPAGGAVK": "P27816[S+80@1073](MAP4)",
        # "AAPEAS[+80]SPPASPLQHLLPGK": "Q96TA1[S+80@691](FAM129B)",
        # "S[+122]DKPDM[+16]AEIEKFDK": "P62328[S+122@2][M+16@7](TMSB4X)",
        # "S[+122]DKPDMAEIEKFDK": "P62328[S+122@2](TMSB4X)",
        # "SLS[+80]LGDKEISR": "Q9UMZ2[S+80@1075](SYNRG)",
        # "DLVQPDKPAS[+80]PK": "Q6PJT7[S+80@515](ZC3H14)",
        # "SPS[+80]PAHLPDDPKVAEK": "Q92615[S+80@601](LARP4B)",
        # "S[+80]IQDLTVTGTEPGQVSSR": "O43318[S+80@439](MAP3K7)",
        # "IHS[+80]PIIR": "O60885[S+80@1117](BRD4)",
        # "TFS[+80]LTEVR": "O95239[S+80@801](KIF4A)",
        # "SLVGS[+80]WLK": "Q6ICG6[S+80@362](KIAA0930)",
        # "S[+80]PPAPGLQPMR": "P15408[S+80@200](FOSL2)",
        # "LAS[+80]PELER": "P17535[S+80@100](JUND)",
        # "IGPLGLS[+80]PK": "P30050[S+80@38](RPL12)",
        # "TPS[+80]IQPSLLPHAAPFAK": "P35658[S+80@1023](NUP214)",
        # "HAS[+80]PILPITEFSDIPR": "P42167[S+80@306](TMPO)",
        # "LIPGPLS[+80]PVAR": "P48634[S+80@1219](PRRC2A)",
        # "LGM[+16]LS[+80]PEGTC[+57]K": "P49327[S+80@207][M+16@205][C+57@212](FASN)",
        # "LGMLS[+80]PEGTC[+57]K": "P49327[S+80@207][C+57@212](FASN)",
        # "ISNLS[+80]PEEEQGLWK": "Q5HYJ3[S+80@193](FAM76B)",
        # "VSMPDVELNLKS[+80]PK": "Q09666[S+80@3426](AHNAK)",
        # "S[+122]DNGELEDKPPAPPVR": "Q13177[S+122@2](PAK2)",
        # "KAYS[+80]FC[+57]GTVEYM[+16]APEVVNR": "Q15418[S+80@221][M+16@229][C+57@223](RPS6KA1)",
        # "KAYS[+80]FC[+57]GTVEYMAPEVVNR": "Q15418[S+80@221][C+57@223](RPS6KA1)",
        # "NDS[+80]WGSFDLR": "Q7Z417[S+80@652](NUFIP2)",
        # "LEVTEIVKPS[+80]PK": "Q7Z6E9[S+80@1179](RBBP6)",
        # "YGS[+80]PPQRDPNWNGER": "O15234[S+80@265](CASC3)",
        # "QDDS[+80]PPRPIIGPALPPGFIK": "Q8IXQ4[S+80@105](GPALPP1)",
        # "SFS[+80]ADNFIGIQR": "Q8N7R7[S+80@344](CCNYL1)",
        # "VLS[+80]PLIIK": "Q8NCN4[S+80@403](RNF169)",
        # "AGS[+80]PDVLR": "Q8NDX6[S+80@44](ZNF740)",
        # "LGPGRPLPTFPTSEC[+57]TS[+80]DVEPDTR": "Q8TDD1[S+80@75][C+57@73](DDX54)",
        # "LAAPSVSHVS[+80]PR": "Q8WXE1[S+80@224](ATRIP)",
        # "VDDDS[+80]LGEFPVTNSR": "Q92785[S+80@142](DPF2)",
        # "NEEPVRS[+80]PERR": "Q92922[S+80@310](SMARCC1)",
        # "LFIIRGS[+80]PQQIDHAK": "Q92945[S+80@480](KHSRP)",
        # "S[+80]IEVENDFLPVEK": "Q96B97[S+80@230](SH3KBP1)",
        # "TAPTLS[+80]PEHWK": "Q96JM3[S+80@405](CHAMP1)",
        # "VLS[+80]PTAAKPSPFEGK": "Q96QC0[S+80@313](PPP1R10)",
        # "SSDQPLTVPVS[+80]PK": "Q9ULW0[S+80@738](TPX2)",
        # "FYETKEESYS[+80]PSKDR": "Q96T23[S+80@473](RSF1)",
        # "SDS[+80]PENKYSDSTGHSK": "Q9BTA9[S+80@64](WAC)",
        # "S[+80]IPLSIK": "Q9C0C9[S+80@515](UBE2O)",
        # "RLS[+80]QSDEDVIR": "Q9H7D7[S+80@121](WDR26)",
        # "ATS[+80]PVKSTTSITDAK": "Q9NQW6[S+80@295](ANLN)",
        # "ALGS[+80]PTKQLLPC[+57]EMAC[+57]NEK": "Q9NR45[S+80@275][C+57@283][C+57@287](NANS)",
        # "YLLGDAPVS[+80]PSSQK": "Q9NYB0[S+80@203](TERF2IP)",
        # "ANS[+80]PEKPPEAGAAHKPR": "Q9UFC0[S+80@212](LRWD1)",
        # "SEVQQPVHPKPLS[+80]PDSR": "Q9UHB6[S+80@362](LIMA1)",
        # "ETPHS[+80]PGVEDAPIAK": "Q9UHB6[S+80@490](LIMA1)",
        # "SQS[+80]PHYFR": "Q9UKJ3[S+80@1035](GPATCH8)",
        # "DRS[+80]SPPPGYIPDELHQVAR": "Q9Y2U5[S+80@163](MAP3K2)",
        # "SPALKS[+80]PLQSVVVR": "Q9Y2W1[S+80@253](THRAP3)",
        # "AFGSGIDIKPGT[+80]PPIAGR": "Q9Y520[T+80@2673](PRRC2C)",
        # "SFS[+80]SQRPVDR": "Q9Y520[S+80@1544](PRRC2C)",
        # "VYT[+80]HEVVTLWYR": "P06493[T+80@161](CDK1)",
        # "SST[+80]PLPTISSSAENTR": "P42167[T+80@160](TMPO)",
        # "QIT[+80]MEELVR": "Q15149[T+80@4030](PLEC)",
        # "TQLWASEPGT[+80]PPLPTSLPSQNPILK": "Q9BXP5[T+80@544](SRRT)",
        # "ALPQT[+80]PRPR": "Q9UQ35[T+80@1492](SRRM2)",
        # "SMS[+80]VDLSHIPLKDPLLFK": "A0JNW5[S+80@935](UHRF1BP1L)",
        # "S[+80]PTGPSNSFLANMGGTVAHK": "Q96I25[S+80@222](RBM17)",
        # "S[+80]LTAHSLLPLAEK": "Q86VI3[S+80@1424](IQGAP3)",
        # "S[+80]FAGNLNTYKR": "Q01813[S+80@386](PFKP)",
        # "HRPS[+80]PPATPPPK": "Q8IYB3[S+80@402](SRRM1)",
        # "LHS[+80]APNLSDLHVVRPK": "O75385[S+80@556](ULK1)",
        # "TLGRRDS[+80]SDDWEIPDGQITVGQR": "P15056[S+80@446](BRAF)",
        # "A[+42]TTATM[+16]ATSGS[+80]AR": "P38919[S+80@12][M+16@7][A+42@2](EIF4A3)",
        # "A[+42]TTATMATSGS[+80]AR": "P38919[S+80@12][A+42@2](EIF4A3)",
        # "IHVSRS[+80]PTRPR": "Q499Z4[S+80@189](ZNF672)",
        # "RPHS[+80]PEKAFSSNPVVR": "Q53F19[S+80@500](NCBP3)",
        # "KPNIFYSGPAS[+80]PARPR": "Q6PL18[S+80@327](ATAD2)",
        # "TEFLDLDNSPLSPPS[+80]PR": "Q8NCF5[S+80@204](NFATC2IP)",
        # "QGSGRES[+80]PSLASR": "Q8WWM7[S+80@339](ATXN2L)",
        # # "TQLWASEPGT[+80]PPLPTSLPSQNPILK": "Q8WWM7[S+80@339](SRRM2)",
        # "LQS[+80]EPESIR": "P09496[S+80@105](CLTA)",
        # "RLIS[+80]PYKK": "O14929[S+80@361](HAT1)",
        # "LLEDS[+80]EESSEETVSR": "O60231[S+80@103](DHX16)",
        # "S[+80]PPAPGLQPM[+16]R": "P15408[S+80@200][M+16@209](FOSL2)",
        # "RRLS[+80]SLR": "P62753[S+80@235](RPS6)",
        # "RLS[+80]ESQLSFRR": "Q96PK6[S+80@618](RBM14)",
        # "RLS[+80]LPGLLSQVSPR": "Q96Q42[S+80@483](ALS2)",
        # "SPDKPGGS[+80]PSASRR": "Q9Y3T9[S+80@56](NOC2L)",
        # "HLPS[+80]PPTLDSIITEYLR": "Q9Y4B6[S+80@1000](DCAF1)",
        # "ST[+80]FHAGQLR": "Q7KZI7[T+80@596](MARK2)",
        # "S[+80]LTNSHLEKK": "Q9H2H9[S+80@52](SLC38A1)",
        # "LQTPNT[+80]FPKR": "Q14978[T+80@610](NOLC1)",
        # "QIT[+80]M[+16]EELVR": "Q15149[M+16@4031][T+80@4030](PLEC)",

        "T[+56]K[+56]QTAR": "P68431[me0K@5](HIST1H3A)",
        "T[+56]K[+70]QTAR": "P68431[meK@5](HIST1H3A)",
        "T[+56]K[+28]QTAR": "P68431[me2K@5](HIST1H3A)",
        "T[+56]K[+42]QTAR-me3K": "P68431[me3K@5](HIST1H3A)",
        "T[+56]K[+42]QTAR-aK": "P68431[aK@5](HIST1H3A)",
        "K[+112.1]STGGK[+56]APR": "P68431[me0K@10][a0K@15](HIST1H3A)",
        "K[+126.1]STGGK[+56]APR": "P68431[meK@10](HIST1H3A)",
        "K[+84.1]STGGK[+56]APR": "P68431[me2K@10](HIST1H3A)",
        "K[+98.1]STGGK[+56]APR": "P68431[me3K@10](HIST1H3A)",
        "K[+98]STGGK[+56]APR": "P68431[aK@10](HIST1H3A)",
        "K[+112.1]STGGK[+42]APR": "P68431[aK@15](HIST1H3A)",
        "K[+126.1]STGGK[+42]APR": "P68431[meK@10][aK@15](HIST1H3A)",
        "K[+84.1]STGGK[+42]APR": "P68431[me2K@10][aK@15](HIST1H3A)",
        "K[+98.1]STGGK[+42]APR": "P68431[me3K@10][aK@15](HIST1H3A)",
        "K[+98]STGGK[+42]APR": "P68431[aK@10][aK@15](HIST1H3A)",
        "K[+112.1]S[+80]TGGK[+56]APR": "P68431[pS@11](HIST1H3A)",
        "K[+126.1]S[+80]TGGK[+56]APR": "P68431[meK@10][pS@11](HIST1H3A)",
        "K[+84.1]S[+80]TGGK[+56]APR": "P68431[me2K@10][pS@11](HIST1H3A)",
        "K[+98.1]S[+80]TGGK[+56]APR": "P68431[me3K@10][pS@11](HIST1H3A)",
        "K[+98]S[+80]TGGK[+56]APR": "P68431[aK@10][pS@11](HIST1H3A)",
        "K[+112.1]S[+80]TGGK[+42]APR": "P68431[pS@11][aK@15](HIST1H3A)",
        "K[+126.1]S[+80]TGGK[+42]APR": "P68431[meK@10][pS@11][aK@15](HIST1H3A)",
        "K[+84.1]S[+80]TGGK[+42]APR": "P68431[me2K@10][pS@11][aK@15](HIST1H3A)",
        "K[+98.1]S[+80]TGGK[+42]APR": "P68431[me3K@10][pS@11][aK@15](HIST1H3A)",
        "K[+98]S[+80]TGGK[+42]APR": "P68431[aK@10][pS@11][aK@15](HIST1H3A)",

        "K[+112.1]QLATK[+56]AAR": "P68431[a0K@19][a0K@24](HIST1H3A)",
        "K[+98]QLATK[+56]AAR": "P68431[aK@19](HIST1H3A)",
        "K[+112.1]QLATK[+42]AAR": "P68431[aK@24](HIST1H3A)",
        "K[+98]QLATK[+42]AAR": "P68431[aK@19][aK@24](HIST1H3A)",

        "K[+226.1]QLATK[+56]AAR": "P68431[ubK@19][a0K@24](HIST1H3A)",
        "K[+112.1]QLATK[+170.1]AAR": "P68431[ubK@24](HIST1H3A)",
        "K[+112.1]SAPATGGVK[+56]K[+56]PHR": "P68431[me0K@28][me0K@37](HIST1H3A)",
        "K[+126.1]SAPATGGVK[+56]K[+56]PHR": "P68431[meK@28][me0K@37](HIST1H3A)",
        "K[+126.1]SAPATGGVK[+70]K[+56]PHR": "P68431[meK@28][meK@37](HIST1H3A)",
        "K[+126.1]SAPATGGVK[+28]K[+56]PHR": "P68431[meK@28][me2K@37](HIST1H3A)",
        "K[+126.1]SAPATGGVK[+42]K[+56]PHR": "P68431[meK@28][me3K@37](HIST1H3A)",
        "K[+84.1]SAPATGGVK[+56]K[+56]PHR": "P68431[me2K@28][me0K@37](HIST1H3A)",
        "K[+84.1]SAPATGGVK[+70]K[+56]PHR": "P68431[me2K@28][meK@37](HIST1H3A)",
        "K[+84.1]SAPATGGVK[+28]K[+56]PHR": "P68431[me2K@28][me2K@37](HIST1H3A)",
        "K[+84.1]SAPATGGVK[+42]K[+56]PHR": "P68431[me2K@28][me3K@37](HIST1H3A)",
        "K[+98.1]SAPATGGVK[+56]K[+56]PHR": "P68431[me3K@28][me0K@37](HIST1H3A)",
        "K[+98.1]SAPATGGVK[+70]K[+56]PHR": "P68431[me3K@28][meK@37](HIST1H3A)",
        "K[+98.1]SAPATGGVK[+28]K[+56]PHR": "P68431[me3K@28][me2K@37](HIST1H3A)",
        "K[+98.1]SAPATGGVK[+42]K[+56]PHR": "P68431[me3K@28][me3K@37](HIST1H3A)",
        "K[+98]SAPATGGVK[+56]K[+56]PHR": "P68431[aK@28][me0K@37](HIST1H3A)",
        "K[+98]SAPATGGVK[+70]K[+56]PHR": "P68431[aK@28][meK@37](HIST1H3A)",
        "K[+98]SAPATGGVK[+28]K[+56]PHR": "P68431[aK@28][me2K@37](HIST1H3A)",
        "K[+98]SAPATGGVK[+42]K[+56]PHR": "P68431[aK@28][me3K@37](HIST1H3A)",
        "K[+112.1]SAPSTGGVK[+56]K[+56]PHR": "P84243[me0K@28][me0K@37](H3F3A)",
        "Y[+56]RPGTVALR": "P68431-NORM(HIST1H3A)",
        "Y[+56]QK[+56]STELLIR": "P68431[me0K@57](HIST1H3A)",
        "E[+56]IAQDFK[+56]TDLR": "P68431[me0K@80](HIST1H3A)",
        "E[+56]IAQDFK[+70]TDLR": "P68431[meK@80](HIST1H3A)",
        "E[+56]IAQDFK[+28]TDLR": "P68431[me2K@80](HIST1H3A)",
        "K[+112.1]SAPATGGVK[+70]K[+56]PHR": "P68431[me0K@28][meK@37](HIST1H3A)",
        "K[+112.1]SAPATGGVK[+28]K[+56]PHR": "P68431[me0K@28][me2K@37](HIST1H3A)",
        "K[+112.1]SAPATGGVK[+42]K[+56]PHR": "P68431[me0K@28][me3K@37](HIST1H3A)",
        "Y[+56]QK[+42]STELLIR": "P68431[aK@57](HIST1H3A)",
        "G[+56]K[+56]GGK[+56]GLGK[+56]GGAK[+56]R": "P62805[a0K@6][a0K@13][a0K@17](HIST1H4A)",
        "G[+56]K[+42]GGK[+56]GLGK[+56]GGAK[+56]R": "P62805[aK@6](HIST1H4A)",
        "G[+56]K[+56]GGK[+56]GLGK[+42]GGAK[+56]R": "P62805[aK@13](HIST1H4A)",
        "G[+56]K[+56]GGK[+56]GLGK[+56]GGAK[+42]R": "P62805[aK@17](HIST1H4A)",
        "G[+56]K[+56]GGK[+42]GLGK[+42]GGAK[+56]R": "P62805[aK@9][aK@13](HIST1H4A)",
        "G[+56]K[+42]GGK[+42]GLGK[+56]GGAK[+56]R": "P62805[aK@6][aK@9](HIST1H4A)",
        "G[+56]K[+42]GGK[+56]GLGK[+56]GGAK[+42]R": "P62805[aK@6][aK@17](HIST1H4A)",
        "G[+56]K[+56]GGK[+56]GLGK[+42]GGAK[+42]R": "P62805[aK@13][aK@17](HIST1H4A)",
        "G[+56]K[+56]GGK[+42]GLGK[+42]GGAK[+42]R": "P62805[aK@9][aK@13][aK@17](HIST1H4A)",
        "G[+56]K[+42]GGK[+42]GLGK[+42]GGAK[+56]R": "P62805[aK@6][aK@9][aK@13](HIST1H4A)",
        "G[+56]K[+42]GGK[+42]GLGK[+56]GGAK[+42]R": "P62805[aK@6][aK@9][aK@17](HIST1H4A)",
        "G[+56]K[+42]GGK[+42]GLGK[+42]GGAK[+42]R": "P62805[aK@6][aK@9][aK@13][aK@17](HIST1H4A)",
        "G[+56]K[+56]GGK[+56]GLGK[+56]GGAK[+70]R": "P62805[meK@17](HIST1H4A)",
        "K[+112.1]VLR": "P62805[me0K@21](HIST1H4A)",
        "K[+126.1]VLR": "P62805[meK@21](HIST1H4A)",
        "K[+84.1]VLR": "P62805[me2K@21](HIST1H4A)",
        "K[+98.1]VLR": "P62805[me3K@21](HIST1H4A)",
        "D[+56]AVTYTEHAK[+56]R": "P62805-NORM(HIST1H4A)",
        "Y[+56]QK[+28]STELLIR": "P68431[me2K@57](HIST1H3A)"
    }
    iter = 0
    for key in peptide2ptm:
        value = peptide2ptm[key]
        print value
        iter += 1

    print iter


    return

@app.route("/api/GCP/downlowdslice/<input_list>")
def GCPSlicingDownload(input_list):
    compounds = []
    results = {}
    for var in input_list.split(","):
        compounds.append(var.upper())

    #
    row_peptide = []
    row_protein = []
    row_gene = []
    row_shorthand = []
    peptide2ptmAndGene = {
        "IYQY[+80]IQSR": "Q13627[Y+80@321](DYRK1A)",
        "TPKDS[+80]PGIPPSANAHQLFR": "P51812[S+80@369](RPS6KA3)",
        "RNS[+80]SEASSGDFLDLK": "Q9UK76[S+80@87](JPT1)",
        "LPLVPES[+80]PRR": "Q86WB0[S+80@321](ZC3HC1)",
        "ANAS[+80]PQKPLDLK": "Q9Y618[S+80@956](NCOR2)",
        "LENS[+80]PLGEALR": "Q9NX40[S+80@108](OCIAD1)",
        "ANS[+80]FVGTAQYVSPELLTEK": "O15530[S+80@241](PDPK1)",
        "TNPPTQKPPS[+80]PPMSGR": "Q8IZP0[S+80@183](ABI1)",
        "SNS[+80]LPHSAVSNAGSK": "Q8TBZ3[S+80@434](WDR20)",
        "VGS[+80]LDNVGHLPAGGAVK": "P27816[S+80@1073](MAP4)",
        "AAPEAS[+80]SPPASPLQHLLPGK": "Q96TA1[S+80@691](FAM129B)",
        "S[+122]DKPDM[+16]AEIEKFDK": "P62328[S+122@2][M+16@7](TMSB4X)",
        "S[+122]DKPDMAEIEKFDK": "P62328[S+122@2](TMSB4X)",
        "SLS[+80]LGDKEISR": "Q9UMZ2[S+80@1075](SYNRG)",
        "DLVQPDKPAS[+80]PK": "Q6PJT7[S+80@515](ZC3H14)",
        "SPS[+80]PAHLPDDPKVAEK": "Q92615[S+80@601](LARP4B)",
        "S[+80]IQDLTVTGTEPGQVSSR": "O43318[S+80@439](MAP3K7)",
        "IHS[+80]PIIR": "O60885[S+80@1117](BRD4)",
        "TFS[+80]LTEVR": "O95239[S+80@801](KIF4A)",
        "SLVGS[+80]WLK": "Q6ICG6[S+80@362](KIAA0930)",
        "S[+80]PPAPGLQPMR": "P15408[S+80@200](FOSL2)",
        "LAS[+80]PELER": "P17535[S+80@100](JUND)",
        "IGPLGLS[+80]PK": "P30050[S+80@38](RPL12)",
        "TPS[+80]IQPSLLPHAAPFAK": "P35658[S+80@1023](NUP214)",
        "HAS[+80]PILPITEFSDIPR": "P42167[S+80@306](TMPO)",
        "LIPGPLS[+80]PVAR": "P48634[S+80@1219](PRRC2A)",
        "LGM[+16]LS[+80]PEGTC[+57]K": "P49327[S+80@207][M+16@205][C+57@212](FASN)",
        "LGMLS[+80]PEGTC[+57]K": "P49327[S+80@207][C+57@212](FASN)",
        "ISNLS[+80]PEEEQGLWK": "Q5HYJ3[S+80@193](FAM76B)",
        "VSMPDVELNLKS[+80]PK": "Q09666[S+80@3426](AHNAK)",
        "S[+122]DNGELEDKPPAPPVR": "Q13177[S+122@2](PAK2)",
        "KAYS[+80]FC[+57]GTVEYM[+16]APEVVNR": "Q15418[S+80@221][M+16@229][C+57@223](RPS6KA1)",
        "KAYS[+80]FC[+57]GTVEYMAPEVVNR": "Q15418[S+80@221][C+57@223](RPS6KA1)",
        "NDS[+80]WGSFDLR": "Q7Z417[S+80@652](NUFIP2)",
        "LEVTEIVKPS[+80]PK": "Q7Z6E9[S+80@1179](RBBP6)",
        "YGS[+80]PPQRDPNWNGER": "O15234[S+80@265](CASC3)",
        "QDDS[+80]PPRPIIGPALPPGFIK": "Q8IXQ4[S+80@105](GPALPP1)",
        "SFS[+80]ADNFIGIQR": "Q8N7R7[S+80@344](CCNYL1)",
        "VLS[+80]PLIIK": "Q8NCN4[S+80@403](RNF169)",
        "AGS[+80]PDVLR": "Q8NDX6[S+80@44](ZNF740)",
        "LGPGRPLPTFPTSEC[+57]TS[+80]DVEPDTR": "Q8TDD1[S+80@75][C+57@73](DDX54)",
        "LAAPSVSHVS[+80]PR": "Q8WXE1[S+80@224](ATRIP)",
        "VDDDS[+80]LGEFPVTNSR": "Q92785[S+80@142](DPF2)",
        "NEEPVRS[+80]PERR": "Q92922[S+80@310](SMARCC1)",
        "LFIIRGS[+80]PQQIDHAK": "Q92945[S+80@480](KHSRP)",
        "S[+80]IEVENDFLPVEK": "Q96B97[S+80@230](SH3KBP1)",
        "TAPTLS[+80]PEHWK": "Q96JM3[S+80@405](CHAMP1)",
        "VLS[+80]PTAAKPSPFEGK": "Q96QC0[S+80@313](PPP1R10)",
        "SSDQPLTVPVS[+80]PK": "Q9ULW0[S+80@738](TPX2)",
        "FYETKEESYS[+80]PSKDR": "Q96T23[S+80@473](RSF1)",
        "SDS[+80]PENKYSDSTGHSK": "Q9BTA9[S+80@64](WAC)",
        "S[+80]IPLSIK": "Q9C0C9[S+80@515](UBE2O)",
        "RLS[+80]QSDEDVIR": "Q9H7D7[S+80@121](WDR26)",
        "ATS[+80]PVKSTTSITDAK": "Q9NQW6[S+80@295](ANLN)",
        "ALGS[+80]PTKQLLPC[+57]EMAC[+57]NEK": "Q9NR45[S+80@275][C+57@283][C+57@287](NANS)",
        "YLLGDAPVS[+80]PSSQK": "Q9NYB0[S+80@203](TERF2IP)",
        "ANS[+80]PEKPPEAGAAHKPR": "Q9UFC0[S+80@212](LRWD1)",
        "SEVQQPVHPKPLS[+80]PDSR": "Q9UHB6[S+80@362](LIMA1)",
        "ETPHS[+80]PGVEDAPIAK": "Q9UHB6[S+80@490](LIMA1)",
        "SQS[+80]PHYFR": "Q9UKJ3[S+80@1035](GPATCH8)",
        "DRS[+80]SPPPGYIPDELHQVAR": "Q9Y2U5[S+80@163](MAP3K2)",
        "SPALKS[+80]PLQSVVVR": "Q9Y2W1[S+80@253](THRAP3)",
        "AFGSGIDIKPGT[+80]PPIAGR": "Q9Y520[T+80@2673](PRRC2C)",
        "SFS[+80]SQRPVDR": "Q9Y520[S+80@1544](PRRC2C)",
        "VYT[+80]HEVVTLWYR": "P06493[T+80@161](CDK1)",
        "SST[+80]PLPTISSSAENTR": "P42167[T+80@160](TMPO)",
        "QIT[+80]MEELVR": "Q15149[T+80@4030](PLEC)",
        "TQLWASEPGT[+80]PPLPTSLPSQNPILK": "Q9BXP5[T+80@544](SRRT)",
        "ALPQT[+80]PRPR": "Q9UQ35[T+80@1492](SRRM2)",
        "SMS[+80]VDLSHIPLKDPLLFK": "A0JNW5[S+80@935](UHRF1BP1L)",
        "S[+80]PTGPSNSFLANMGGTVAHK": "Q96I25[S+80@222](RBM17)",
        "S[+80]LTAHSLLPLAEK": "Q86VI3[S+80@1424](IQGAP3)",
        "S[+80]FAGNLNTYKR": "Q01813[S+80@386](PFKP)",
        "HRPS[+80]PPATPPPK": "Q8IYB3[S+80@402](SRRM1)",
        "LHS[+80]APNLSDLHVVRPK": "O75385[S+80@556](ULK1)",
        "TLGRRDS[+80]SDDWEIPDGQITVGQR": "P15056[S+80@446](BRAF)",
        "A[+42]TTATM[+16]ATSGS[+80]AR": "P38919[S+80@12][M+16@7][A+42@2](EIF4A3)",
        "A[+42]TTATMATSGS[+80]AR": "P38919[S+80@12][A+42@2](EIF4A3)",
        "IHVSRS[+80]PTRPR": "Q499Z4[S+80@189](ZNF672)",
        "RPHS[+80]PEKAFSSNPVVR": "Q53F19[S+80@500](NCBP3)",
        "KPNIFYSGPAS[+80]PARPR": "Q6PL18[S+80@327](ATAD2)",
        "TEFLDLDNSPLSPPS[+80]PR": "Q8NCF5[S+80@204](NFATC2IP)",
        "QGSGRES[+80]PSLASR": "Q8WWM7[S+80@339](ATXN2L)",
        # "TQLWASEPGT[+80]PPLPTSLPSQNPILK": "Q8WWM7[S+80@339](SRRM2)",
        "LQS[+80]EPESIR": "P09496[S+80@105](CLTA)",
        "RLIS[+80]PYKK": "O14929[S+80@361](HAT1)",
        "LLEDS[+80]EESSEETVSR": "O60231[S+80@103](DHX16)",
        "S[+80]PPAPGLQPM[+16]R": "P15408[S+80@200][M+16@209](FOSL2)",
        "RRLS[+80]SLR": "P62753[S+80@235](RPS6)",
        "RLS[+80]ESQLSFRR": "Q96PK6[S+80@618](RBM14)",
        "RLS[+80]LPGLLSQVSPR": "Q96Q42[S+80@483](ALS2)",
        "SPDKPGGS[+80]PSASRR": "Q9Y3T9[S+80@56](NOC2L)",
        "HLPS[+80]PPTLDSIITEYLR": "Q9Y4B6[S+80@1000](DCAF1)",
        "ST[+80]FHAGQLR": "Q7KZI7[T+80@596](MARK2)",
        "S[+80]LTNSHLEKK": "Q9H2H9[S+80@52](SLC38A1)",
        "LQTPNT[+80]FPKR": "Q14978[T+80@610](NOLC1)",
        "QIT[+80]M[+16]EELVR": "Q15149[M+16@4031][T+80@4030](PLEC)",

        "T[+56]K[+56]QTAR": "P68431[me0K@5](HIST1H3A)",
        "T[+56]K[+70]QTAR": "P68431[meK@5](HIST1H3A)",
        "T[+56]K[+28]QTAR": "P68431[me2K@5](HIST1H3A)",
        "T[+56]K[+42]QTAR-me3K": "P68431[me3K@5](HIST1H3A)",
        "T[+56]K[+42]QTAR-aK": "P68431[aK@5](HIST1H3A)",
        "K[+112.1]STGGK[+56]APR": "P68431[me0K@10][a0K@15](HIST1H3A)",
        "K[+126.1]STGGK[+56]APR": "P68431[meK@10](HIST1H3A)",
        "K[+84.1]STGGK[+56]APR": "P68431[me2K@10](HIST1H3A)",
        "K[+98.1]STGGK[+56]APR": "P68431[me3K@10](HIST1H3A)",
        "K[+98]STGGK[+56]APR": "P68431[aK@10](HIST1H3A)",
        "K[+112.1]STGGK[+42]APR": "P68431[aK@15](HIST1H3A)",
        "K[+126.1]STGGK[+42]APR": "P68431[meK@10][aK@15](HIST1H3A)",
        "K[+84.1]STGGK[+42]APR": "P68431[me2K@10][aK@15](HIST1H3A)",
        "K[+98.1]STGGK[+42]APR": "P68431[me3K@10][aK@15](HIST1H3A)",
        "K[+98]STGGK[+42]APR": "P68431[aK@10][aK@15](HIST1H3A)",
        "K[+112.1]S[+80]TGGK[+56]APR": "P68431[pS@11](HIST1H3A)",
        "K[+126.1]S[+80]TGGK[+56]APR": "P68431[meK@10][pS@11](HIST1H3A)",
        "K[+84.1]S[+80]TGGK[+56]APR": "P68431[me2K@10][pS@11](HIST1H3A)",
        "K[+98.1]S[+80]TGGK[+56]APR": "P68431[me3K@10][pS@11](HIST1H3A)",
        "K[+98]S[+80]TGGK[+56]APR": "P68431[aK@10][pS@11](HIST1H3A)",
        "K[+112.1]S[+80]TGGK[+42]APR": "P68431[pS@11][aK@15](HIST1H3A)",
        "K[+126.1]S[+80]TGGK[+42]APR": "P68431[meK@10][pS@11][aK@15](HIST1H3A)",
        "K[+84.1]S[+80]TGGK[+42]APR": "P68431[me2K@10][pS@11][aK@15](HIST1H3A)",
        "K[+98.1]S[+80]TGGK[+42]APR": "P68431[me3K@10][pS@11][aK@15](HIST1H3A)",
        "K[+98]S[+80]TGGK[+42]APR": "P68431[aK@10][pS@11][aK@15](HIST1H3A)",

        "K[+112.1]QLATK[+56]AAR": "P68431[a0K@19][a0K@24](HIST1H3A)",
        "K[+98]QLATK[+56]AAR": "P68431[aK@19](HIST1H3A)",
        "K[+112.1]QLATK[+42]AAR": "P68431[aK@24](HIST1H3A)",
        "K[+98]QLATK[+42]AAR": "P68431[aK@19][aK@24](HIST1H3A)",

        "K[+226.1]QLATK[+56]AAR": "P68431[ubK@19][a0K@24](HIST1H3A)",
        "K[+112.1]QLATK[+170.1]AAR": "P68431[ubK@24](HIST1H3A)",
        "K[+112.1]SAPATGGVK[+56]K[+56]PHR": "P68431[me0K@28][me0K@37](HIST1H3A)",
        "K[+126.1]SAPATGGVK[+56]K[+56]PHR": "P68431[meK@28][me0K@37](HIST1H3A)",
        "K[+126.1]SAPATGGVK[+70]K[+56]PHR": "P68431[meK@28][meK@37](HIST1H3A)",
        "K[+126.1]SAPATGGVK[+28]K[+56]PHR": "P68431[meK@28][me2K@37](HIST1H3A)",
        "K[+126.1]SAPATGGVK[+42]K[+56]PHR": "P68431[meK@28][me3K@37](HIST1H3A)",
        "K[+84.1]SAPATGGVK[+56]K[+56]PHR": "P68431[me2K@28][me0K@37](HIST1H3A)",
        "K[+84.1]SAPATGGVK[+70]K[+56]PHR": "P68431[me2K@28][meK@37](HIST1H3A)",
        "K[+84.1]SAPATGGVK[+28]K[+56]PHR": "P68431[me2K@28][me2K@37](HIST1H3A)",
        "K[+84.1]SAPATGGVK[+42]K[+56]PHR": "P68431[me2K@28][me3K@37](HIST1H3A)",
        "K[+98.1]SAPATGGVK[+56]K[+56]PHR": "P68431[me3K@28][me0K@37](HIST1H3A)",
        "K[+98.1]SAPATGGVK[+70]K[+56]PHR": "P68431[me3K@28][meK@37](HIST1H3A)",
        "K[+98.1]SAPATGGVK[+28]K[+56]PHR": "P68431[me3K@28][me2K@37](HIST1H3A)",
        "K[+98.1]SAPATGGVK[+42]K[+56]PHR": "P68431[me3K@28][me3K@37](HIST1H3A)",
        "K[+98]SAPATGGVK[+56]K[+56]PHR": "P68431[aK@28][me0K@37](HIST1H3A)",
        "K[+98]SAPATGGVK[+70]K[+56]PHR": "P68431[aK@28][meK@37](HIST1H3A)",
        "K[+98]SAPATGGVK[+28]K[+56]PHR": "P68431[aK@28][me2K@37](HIST1H3A)",
        "K[+98]SAPATGGVK[+42]K[+56]PHR": "P68431[aK@28][me3K@37](HIST1H3A)",
        "K[+112.1]SAPSTGGVK[+56]K[+56]PHR": "P84243[me0K@28][me0K@37](H3F3A)",
        "Y[+56]RPGTVALR": "P68431-NORM(HIST1H3A)",
        "Y[+56]QK[+56]STELLIR": "P68431[me0K@57](HIST1H3A)",
        "E[+56]IAQDFK[+56]TDLR": "P68431[me0K@80](HIST1H3A)",
        "E[+56]IAQDFK[+70]TDLR": "P68431[meK@80](HIST1H3A)",
        "E[+56]IAQDFK[+28]TDLR": "P68431[me2K@80](HIST1H3A)",
        "K[+112.1]SAPATGGVK[+70]K[+56]PHR": "P68431[me0K@28][meK@37](HIST1H3A)",
        "K[+112.1]SAPATGGVK[+28]K[+56]PHR": "P68431[me0K@28][me2K@37](HIST1H3A)",
        "K[+112.1]SAPATGGVK[+42]K[+56]PHR": "P68431[me0K@28][me3K@37](HIST1H3A)",
        "Y[+56]QK[+42]STELLIR": "P68431[aK@57](HIST1H3A)",
        "G[+56]K[+56]GGK[+56]GLGK[+56]GGAK[+56]R": "P62805[a0K@6][a0K@13][a0K@17](HIST1H4A)",
        "G[+56]K[+42]GGK[+56]GLGK[+56]GGAK[+56]R": "P62805[aK@6](HIST1H4A)",
        "G[+56]K[+56]GGK[+56]GLGK[+42]GGAK[+56]R": "P62805[aK@13](HIST1H4A)",
        "G[+56]K[+56]GGK[+56]GLGK[+56]GGAK[+42]R": "P62805[aK@17](HIST1H4A)",
        "G[+56]K[+56]GGK[+42]GLGK[+42]GGAK[+56]R": "P62805[aK@9][aK@13](HIST1H4A)",
        "G[+56]K[+42]GGK[+42]GLGK[+56]GGAK[+56]R": "P62805[aK@6][aK@9](HIST1H4A)",
        "G[+56]K[+42]GGK[+56]GLGK[+56]GGAK[+42]R": "P62805[aK@6][aK@17](HIST1H4A)",
        "G[+56]K[+56]GGK[+56]GLGK[+42]GGAK[+42]R": "P62805[aK@13][aK@17](HIST1H4A)",
        "G[+56]K[+56]GGK[+42]GLGK[+42]GGAK[+42]R": "P62805[aK@9][aK@13][aK@17](HIST1H4A)",
        "G[+56]K[+42]GGK[+42]GLGK[+42]GGAK[+56]R": "P62805[aK@6][aK@9][aK@13](HIST1H4A)",
        "G[+56]K[+42]GGK[+42]GLGK[+56]GGAK[+42]R": "P62805[aK@6][aK@9][aK@17](HIST1H4A)",
        "G[+56]K[+42]GGK[+42]GLGK[+42]GGAK[+42]R": "P62805[aK@6][aK@9][aK@13][aK@17](HIST1H4A)",
        "G[+56]K[+56]GGK[+56]GLGK[+56]GGAK[+70]R": "P62805[meK@17](HIST1H4A)",
        "K[+112.1]VLR": "P62805[me0K@21](HIST1H4A)",
        "K[+126.1]VLR": "P62805[meK@21](HIST1H4A)",
        "K[+84.1]VLR": "P62805[me2K@21](HIST1H4A)",
        "K[+98.1]VLR": "P62805[me3K@21](HIST1H4A)",
        "D[+56]AVTYTEHAK[+56]R": "P62805-NORM(HIST1H4A)",
        "Y[+56]QK[+28]STELLIR": "P68431[me2K@57](HIST1H3A)"
    }
    meta_data = pd.read_csv('static/data/GCP-all-plates-metaData.csv', sep=',', header=None)
    #print meta_data

    for i in range(0, 79):

        peptide = str(meta_data[4][28 + i])
        if(peptide == 'T[+56]K[+42]QTAR'):
            print str(meta_data[2][28 + i])
            if(str(meta_data[2][28 + i]) == 'BI10006'):

                peptide = "T[+56]K[+42]QTAR-me3K"
            else:
                peptide = "T[+56]K[+42]QTAR-aK"
        print peptide
        row_peptide.append(peptide)

        shorthand = peptide2ptmAndGene[peptide]

        row_shorthand.append('PTM Proteins: ' + shorthand)
        row_protein.append('Protein: ' + str(meta_data[9][28 + i]))
        row_gene.append('Gene: ' + str(meta_data[6][28 + i]))



    col_cell_line = []
    col_pert_dose = []
    col_pert_time = []
    col_pert_name = []
    col_pert_signature = []

    col_pert_type = []
    first_loop = True




    # for i in range(0, 96):
    #     row_peptide.append('Peptide: ' + str(P100_data[6][29 + i]))
    #     shorthand = petide2ptmAndGene[str(P100_data[6][29 + i])]
    #     row_shorthand.append('PTM Proteins: ' + shorthand)
    #     row_protein.append('Protein: ' + str(P100_data[11][29 + i]))
    #     row_gene.append('Gene: ' + str(P100_data[2][29 + i]))

    # open('static/data/P100_processed_perturb_for_clustergramm_slicing.json') as f:



    with open('static/data/GCP_processed_perturb_for_clustergramm_slicing.json') as f:
        data = json.load(f)
        cp_iter = 1
        for var in compounds:
            cp = str(var).upper()
            print cp

            if cp in data:
                for iter in data[cp]:

                    print cp_iter
                    cp_iter += 1
                    col_cell_line.append('Cell Line: ' + iter["Cell_line"])
                    col_pert_dose.append('Dose: ' + iter["Dose"])
                    col_pert_time.append('Time: ' + iter["Time"])
                    col_pert_name.append(iter["Perturbations"]+"/Dose/"+iter["Dose"]+"/Time/"+iter["Time"]+"/CellLine/"+iter["Cell_line"]+"/Signature ID/"+iter["Sig_id"])
                    col_pert_signature.append('Signature ID: ' + iter["Sig_id"])
                    col_pert_type.append('Type: ' + iter["Type"])

                    x1 = np.array(iter["data"])

                    # print "end of loop1"
                    x2 = x1[:, np.newaxis]
                    # print "end of loop2"
                    my_data_final = np.asarray(x2)
                    # print "end of loop3"
                    if first_loop:
                        my_data_con = my_data_final
                        first_loop = False

                    else:
                        my_data_con = np.hstack((my_data_con, my_data_final))


        # for i in range(0, 96):
        #     row_peptide.append('Peptide: ' + str(P100_data[6][29 + i]))
        #     shorthand = petide2ptmAndGene[str(P100_data[6][29 + i])]
        #     row_shorthand.append('PTM Proteins: ' + shorthand)
        #     row_protein.append('Protein: ' + str(P100_data[11][29 + i]))
        #     row_gene.append('Gene: ' + str(P100_data[2][29 + i]))

        arrays_rows2 = [np.array(row_peptide)]
        arrays_columns2 = [np.array(col_pert_name)]
        # print arrays_columns2
        # arrays_rows2 = np.load("static/data/P100_rowMetaDat")

        tuples_rows2 = list(zip(*arrays_rows2))
        # print tuples_rows2
        tuples_columns2 = list(zip(*arrays_columns2))
        # print tuples_columns2
        rows_labels2 = pd.MultiIndex.from_tuples(tuples_rows2)
        columns_labels2 = pd.MultiIndex.from_tuples(tuples_columns2)

        # my_data_con = my_data_con.fillna(0)
        print "after na to zero"
        a2 = np.nan_to_num(my_data_con)
        # a2 = my_data_con
        aa = a2.astype(np.float)

        print "a2.shape"
        print a2.shape
        print "columns_labels2.shape"
        print columns_labels2.shape
        print "rows_labels2.shape"
        print rows_labels2.shape
        # {'nop': row1, 'o0p': row2, 'zaz': row3, 'zax': row4, 'oof': row5, 'oye': row6}
        # df = pd.DataFrame(a, index=rows_labels, columns=columns_labels)
        df2 = pd.DataFrame(aa, index=rows_labels2, columns=columns_labels2)
        jsonRet = df2.to_dict()

        print jsonRet
        return df2.to_json(orient='split')



            #print results
    # print "finished"
    # GCP_all_json = loadFile.make_json_from_txt(df2)

@app.route("/api/P100/downlowdslice/<input_list>")
def P100SlicingDownload(input_list):
    compounds = []
    results = {}
    for var in input_list.split(","):
        compounds.append(var.upper())

    #
    row_peptide = []
    row_protein = []
    row_gene = []
    row_shorthand = []
    peptide2ptmAndGene = {
        "IYQY[+80]IQSR": "Q13627[Y+80@321](DYRK1A)",
        "TPKDS[+80]PGIPPSANAHQLFR": "P51812[S+80@369](RPS6KA3)",
        "RNS[+80]SEASSGDFLDLK": "Q9UK76[S+80@87](JPT1)",
        "LPLVPES[+80]PRR": "Q86WB0[S+80@321](ZC3HC1)",
        "ANAS[+80]PQKPLDLK": "Q9Y618[S+80@956](NCOR2)",
        "LENS[+80]PLGEALR": "Q9NX40[S+80@108](OCIAD1)",
        "ANS[+80]FVGTAQYVSPELLTEK": "O15530[S+80@241](PDPK1)",
        "TNPPTQKPPS[+80]PPMSGR": "Q8IZP0[S+80@183](ABI1)",
        "SNS[+80]LPHSAVSNAGSK": "Q8TBZ3[S+80@434](WDR20)",
        "VGS[+80]LDNVGHLPAGGAVK": "P27816[S+80@1073](MAP4)",
        "AAPEAS[+80]SPPASPLQHLLPGK": "Q96TA1[S+80@691](FAM129B)",
        "S[+122]DKPDM[+16]AEIEKFDK": "P62328[S+122@2][M+16@7](TMSB4X)",
        "S[+122]DKPDMAEIEKFDK": "P62328[S+122@2](TMSB4X)",
        "SLS[+80]LGDKEISR": "Q9UMZ2[S+80@1075](SYNRG)",
        "DLVQPDKPAS[+80]PK": "Q6PJT7[S+80@515](ZC3H14)",
        "SPS[+80]PAHLPDDPKVAEK": "Q92615[S+80@601](LARP4B)",
        "S[+80]IQDLTVTGTEPGQVSSR": "O43318[S+80@439](MAP3K7)",
        "IHS[+80]PIIR": "O60885[S+80@1117](BRD4)",
        "TFS[+80]LTEVR": "O95239[S+80@801](KIF4A)",
        "SLVGS[+80]WLK": "Q6ICG6[S+80@362](KIAA0930)",
        "S[+80]PPAPGLQPMR": "P15408[S+80@200](FOSL2)",
        "LAS[+80]PELER": "P17535[S+80@100](JUND)",
        "IGPLGLS[+80]PK": "P30050[S+80@38](RPL12)",
        "TPS[+80]IQPSLLPHAAPFAK": "P35658[S+80@1023](NUP214)",
        "HAS[+80]PILPITEFSDIPR": "P42167[S+80@306](TMPO)",
        "LIPGPLS[+80]PVAR": "P48634[S+80@1219](PRRC2A)",
        "LGM[+16]LS[+80]PEGTC[+57]K": "P49327[S+80@207][M+16@205][C+57@212](FASN)",
        "LGMLS[+80]PEGTC[+57]K": "P49327[S+80@207][C+57@212](FASN)",
        "ISNLS[+80]PEEEQGLWK": "Q5HYJ3[S+80@193](FAM76B)",
        "VSMPDVELNLKS[+80]PK": "Q09666[S+80@3426](AHNAK)",
        "S[+122]DNGELEDKPPAPPVR": "Q13177[S+122@2](PAK2)",
        "KAYS[+80]FC[+57]GTVEYM[+16]APEVVNR": "Q15418[S+80@221][M+16@229][C+57@223](RPS6KA1)",
        "KAYS[+80]FC[+57]GTVEYMAPEVVNR": "Q15418[S+80@221][C+57@223](RPS6KA1)",
        "NDS[+80]WGSFDLR": "Q7Z417[S+80@652](NUFIP2)",
        "LEVTEIVKPS[+80]PK": "Q7Z6E9[S+80@1179](RBBP6)",
        "YGS[+80]PPQRDPNWNGER": "O15234[S+80@265](CASC3)",
        "QDDS[+80]PPRPIIGPALPPGFIK": "Q8IXQ4[S+80@105](GPALPP1)",
        "SFS[+80]ADNFIGIQR": "Q8N7R7[S+80@344](CCNYL1)",
        "VLS[+80]PLIIK": "Q8NCN4[S+80@403](RNF169)",
        "AGS[+80]PDVLR": "Q8NDX6[S+80@44](ZNF740)",
        "LGPGRPLPTFPTSEC[+57]TS[+80]DVEPDTR": "Q8TDD1[S+80@75][C+57@73](DDX54)",
        "LAAPSVSHVS[+80]PR": "Q8WXE1[S+80@224](ATRIP)",
        "VDDDS[+80]LGEFPVTNSR": "Q92785[S+80@142](DPF2)",
        "NEEPVRS[+80]PERR": "Q92922[S+80@310](SMARCC1)",
        "LFIIRGS[+80]PQQIDHAK": "Q92945[S+80@480](KHSRP)",
        "S[+80]IEVENDFLPVEK": "Q96B97[S+80@230](SH3KBP1)",
        "TAPTLS[+80]PEHWK": "Q96JM3[S+80@405](CHAMP1)",
        "VLS[+80]PTAAKPSPFEGK": "Q96QC0[S+80@313](PPP1R10)",
        "SSDQPLTVPVS[+80]PK": "Q9ULW0[S+80@738](TPX2)",
        "FYETKEESYS[+80]PSKDR": "Q96T23[S+80@473](RSF1)",
        "SDS[+80]PENKYSDSTGHSK": "Q9BTA9[S+80@64](WAC)",
        "S[+80]IPLSIK": "Q9C0C9[S+80@515](UBE2O)",
        "RLS[+80]QSDEDVIR": "Q9H7D7[S+80@121](WDR26)",
        "ATS[+80]PVKSTTSITDAK": "Q9NQW6[S+80@295](ANLN)",
        "ALGS[+80]PTKQLLPC[+57]EMAC[+57]NEK": "Q9NR45[S+80@275][C+57@283][C+57@287](NANS)",
        "YLLGDAPVS[+80]PSSQK": "Q9NYB0[S+80@203](TERF2IP)",
        "ANS[+80]PEKPPEAGAAHKPR": "Q9UFC0[S+80@212](LRWD1)",
        "SEVQQPVHPKPLS[+80]PDSR": "Q9UHB6[S+80@362](LIMA1)",
        "ETPHS[+80]PGVEDAPIAK": "Q9UHB6[S+80@490](LIMA1)",
        "SQS[+80]PHYFR": "Q9UKJ3[S+80@1035](GPATCH8)",
        "DRS[+80]SPPPGYIPDELHQVAR": "Q9Y2U5[S+80@163](MAP3K2)",
        "SPALKS[+80]PLQSVVVR": "Q9Y2W1[S+80@253](THRAP3)",
        "AFGSGIDIKPGT[+80]PPIAGR": "Q9Y520[T+80@2673](PRRC2C)",
        "SFS[+80]SQRPVDR": "Q9Y520[S+80@1544](PRRC2C)",
        "VYT[+80]HEVVTLWYR": "P06493[T+80@161](CDK1)",
        "SST[+80]PLPTISSSAENTR": "P42167[T+80@160](TMPO)",
        "QIT[+80]MEELVR": "Q15149[T+80@4030](PLEC)",
        "TQLWASEPGT[+80]PPLPTSLPSQNPILK": "Q9BXP5[T+80@544](SRRT)",
        "ALPQT[+80]PRPR": "Q9UQ35[T+80@1492](SRRM2)",
        "SMS[+80]VDLSHIPLKDPLLFK": "A0JNW5[S+80@935](UHRF1BP1L)",
        "S[+80]PTGPSNSFLANMGGTVAHK": "Q96I25[S+80@222](RBM17)",
        "S[+80]LTAHSLLPLAEK": "Q86VI3[S+80@1424](IQGAP3)",
        "S[+80]FAGNLNTYKR": "Q01813[S+80@386](PFKP)",
        "HRPS[+80]PPATPPPK": "Q8IYB3[S+80@402](SRRM1)",
        "LHS[+80]APNLSDLHVVRPK": "O75385[S+80@556](ULK1)",
        "TLGRRDS[+80]SDDWEIPDGQITVGQR": "P15056[S+80@446](BRAF)",
        "A[+42]TTATM[+16]ATSGS[+80]AR": "P38919[S+80@12][M+16@7][A+42@2](EIF4A3)",
        "A[+42]TTATMATSGS[+80]AR": "P38919[S+80@12][A+42@2](EIF4A3)",
        "IHVSRS[+80]PTRPR": "Q499Z4[S+80@189](ZNF672)",
        "RPHS[+80]PEKAFSSNPVVR": "Q53F19[S+80@500](NCBP3)",
        "KPNIFYSGPAS[+80]PARPR": "Q6PL18[S+80@327](ATAD2)",
        "TEFLDLDNSPLSPPS[+80]PR": "Q8NCF5[S+80@204](NFATC2IP)",
        "QGSGRES[+80]PSLASR": "Q8WWM7[S+80@339](ATXN2L)",
        # "TQLWASEPGT[+80]PPLPTSLPSQNPILK": "Q8WWM7[S+80@339](SRRM2)",
        "LQS[+80]EPESIR": "P09496[S+80@105](CLTA)",
        "RLIS[+80]PYKK": "O14929[S+80@361](HAT1)",
        "LLEDS[+80]EESSEETVSR": "O60231[S+80@103](DHX16)",
        "S[+80]PPAPGLQPM[+16]R": "P15408[S+80@200][M+16@209](FOSL2)",
        "RRLS[+80]SLR": "P62753[S+80@235](RPS6)",
        "RLS[+80]ESQLSFRR": "Q96PK6[S+80@618](RBM14)",
        "RLS[+80]LPGLLSQVSPR": "Q96Q42[S+80@483](ALS2)",
        "SPDKPGGS[+80]PSASRR": "Q9Y3T9[S+80@56](NOC2L)",
        "HLPS[+80]PPTLDSIITEYLR": "Q9Y4B6[S+80@1000](DCAF1)",
        "ST[+80]FHAGQLR": "Q7KZI7[T+80@596](MARK2)",
        "S[+80]LTNSHLEKK": "Q9H2H9[S+80@52](SLC38A1)",
        "LQTPNT[+80]FPKR": "Q14978[T+80@610](NOLC1)",
        "QIT[+80]M[+16]EELVR": "Q15149[M+16@4031][T+80@4030](PLEC)",

        "T[+56]K[+56]QTAR": "P68431[me0K@5](HIST1H3A)",
        "T[+56]K[+70]QTAR": "P68431[meK@5](HIST1H3A)",
        "T[+56]K[+28]QTAR": "P68431[me2K@5](HIST1H3A)",
        "T[+56]K[+42]QTAR-me3K": "P68431[me3K@5](HIST1H3A)",
        "T[+56]K[+42]QTAR-aK": "P68431[aK@5](HIST1H3A)",
        "K[+112.1]STGGK[+56]APR": "P68431[me0K@10][a0K@15](HIST1H3A)",
        "K[+126.1]STGGK[+56]APR": "P68431[meK@10](HIST1H3A)",
        "K[+84.1]STGGK[+56]APR": "P68431[me2K@10](HIST1H3A)",
        "K[+98.1]STGGK[+56]APR": "P68431[me3K@10](HIST1H3A)",
        "K[+98]STGGK[+56]APR": "P68431[aK@10](HIST1H3A)",
        "K[+112.1]STGGK[+42]APR": "P68431[aK@15](HIST1H3A)",
        "K[+126.1]STGGK[+42]APR": "P68431[meK@10][aK@15](HIST1H3A)",
        "K[+84.1]STGGK[+42]APR": "P68431[me2K@10][aK@15](HIST1H3A)",
        "K[+98.1]STGGK[+42]APR": "P68431[me3K@10][aK@15](HIST1H3A)",
        "K[+98]STGGK[+42]APR": "P68431[aK@10][aK@15](HIST1H3A)",
        "K[+112.1]S[+80]TGGK[+56]APR": "P68431[pS@11](HIST1H3A)",
        "K[+126.1]S[+80]TGGK[+56]APR": "P68431[meK@10][pS@11](HIST1H3A)",
        "K[+84.1]S[+80]TGGK[+56]APR": "P68431[me2K@10][pS@11](HIST1H3A)",
        "K[+98.1]S[+80]TGGK[+56]APR": "P68431[me3K@10][pS@11](HIST1H3A)",
        "K[+98]S[+80]TGGK[+56]APR": "P68431[aK@10][pS@11](HIST1H3A)",
        "K[+112.1]S[+80]TGGK[+42]APR": "P68431[pS@11][aK@15](HIST1H3A)",
        "K[+126.1]S[+80]TGGK[+42]APR": "P68431[meK@10][pS@11][aK@15](HIST1H3A)",
        "K[+84.1]S[+80]TGGK[+42]APR": "P68431[me2K@10][pS@11][aK@15](HIST1H3A)",
        "K[+98.1]S[+80]TGGK[+42]APR": "P68431[me3K@10][pS@11][aK@15](HIST1H3A)",
        "K[+98]S[+80]TGGK[+42]APR": "P68431[aK@10][pS@11][aK@15](HIST1H3A)",

        "K[+112.1]QLATK[+56]AAR": "P68431[a0K@19][a0K@24](HIST1H3A)",
        "K[+98]QLATK[+56]AAR": "P68431[aK@19](HIST1H3A)",
        "K[+112.1]QLATK[+42]AAR": "P68431[aK@24](HIST1H3A)",
        "K[+98]QLATK[+42]AAR": "P68431[aK@19][aK@24](HIST1H3A)",

        "K[+226.1]QLATK[+56]AAR": "P68431[ubK@19][a0K@24](HIST1H3A)",
        "K[+112.1]QLATK[+170.1]AAR": "P68431[ubK@24](HIST1H3A)",
        "K[+112.1]SAPATGGVK[+56]K[+56]PHR": "P68431[me0K@28][me0K@37](HIST1H3A)",
        "K[+126.1]SAPATGGVK[+56]K[+56]PHR": "P68431[meK@28][me0K@37](HIST1H3A)",
        "K[+126.1]SAPATGGVK[+70]K[+56]PHR": "P68431[meK@28][meK@37](HIST1H3A)",
        "K[+126.1]SAPATGGVK[+28]K[+56]PHR": "P68431[meK@28][me2K@37](HIST1H3A)",
        "K[+126.1]SAPATGGVK[+42]K[+56]PHR": "P68431[meK@28][me3K@37](HIST1H3A)",
        "K[+84.1]SAPATGGVK[+56]K[+56]PHR": "P68431[me2K@28][me0K@37](HIST1H3A)",
        "K[+84.1]SAPATGGVK[+70]K[+56]PHR": "P68431[me2K@28][meK@37](HIST1H3A)",
        "K[+84.1]SAPATGGVK[+28]K[+56]PHR": "P68431[me2K@28][me2K@37](HIST1H3A)",
        "K[+84.1]SAPATGGVK[+42]K[+56]PHR": "P68431[me2K@28][me3K@37](HIST1H3A)",
        "K[+98.1]SAPATGGVK[+56]K[+56]PHR": "P68431[me3K@28][me0K@37](HIST1H3A)",
        "K[+98.1]SAPATGGVK[+70]K[+56]PHR": "P68431[me3K@28][meK@37](HIST1H3A)",
        "K[+98.1]SAPATGGVK[+28]K[+56]PHR": "P68431[me3K@28][me2K@37](HIST1H3A)",
        "K[+98.1]SAPATGGVK[+42]K[+56]PHR": "P68431[me3K@28][me3K@37](HIST1H3A)",
        "K[+98]SAPATGGVK[+56]K[+56]PHR": "P68431[aK@28][me0K@37](HIST1H3A)",
        "K[+98]SAPATGGVK[+70]K[+56]PHR": "P68431[aK@28][meK@37](HIST1H3A)",
        "K[+98]SAPATGGVK[+28]K[+56]PHR": "P68431[aK@28][me2K@37](HIST1H3A)",
        "K[+98]SAPATGGVK[+42]K[+56]PHR": "P68431[aK@28][me3K@37](HIST1H3A)",
        "K[+112.1]SAPSTGGVK[+56]K[+56]PHR": "P84243[me0K@28][me0K@37](H3F3A)",
        "Y[+56]RPGTVALR": "P68431-NORM(HIST1H3A)",
        "Y[+56]QK[+56]STELLIR": "P68431[me0K@57](HIST1H3A)",
        "E[+56]IAQDFK[+56]TDLR": "P68431[me0K@80](HIST1H3A)",
        "E[+56]IAQDFK[+70]TDLR": "P68431[meK@80](HIST1H3A)",
        "E[+56]IAQDFK[+28]TDLR": "P68431[me2K@80](HIST1H3A)",
        "K[+112.1]SAPATGGVK[+70]K[+56]PHR": "P68431[me0K@28][meK@37](HIST1H3A)",
        "K[+112.1]SAPATGGVK[+28]K[+56]PHR": "P68431[me0K@28][me2K@37](HIST1H3A)",
        "K[+112.1]SAPATGGVK[+42]K[+56]PHR": "P68431[me0K@28][me3K@37](HIST1H3A)",
        "Y[+56]QK[+42]STELLIR": "P68431[aK@57](HIST1H3A)",
        "G[+56]K[+56]GGK[+56]GLGK[+56]GGAK[+56]R": "P62805[a0K@6][a0K@13][a0K@17](HIST1H4A)",
        "G[+56]K[+42]GGK[+56]GLGK[+56]GGAK[+56]R": "P62805[aK@6](HIST1H4A)",
        "G[+56]K[+56]GGK[+56]GLGK[+42]GGAK[+56]R": "P62805[aK@13](HIST1H4A)",
        "G[+56]K[+56]GGK[+56]GLGK[+56]GGAK[+42]R": "P62805[aK@17](HIST1H4A)",
        "G[+56]K[+56]GGK[+42]GLGK[+42]GGAK[+56]R": "P62805[aK@9][aK@13](HIST1H4A)",
        "G[+56]K[+42]GGK[+42]GLGK[+56]GGAK[+56]R": "P62805[aK@6][aK@9](HIST1H4A)",
        "G[+56]K[+42]GGK[+56]GLGK[+56]GGAK[+42]R": "P62805[aK@6][aK@17](HIST1H4A)",
        "G[+56]K[+56]GGK[+56]GLGK[+42]GGAK[+42]R": "P62805[aK@13][aK@17](HIST1H4A)",
        "G[+56]K[+56]GGK[+42]GLGK[+42]GGAK[+42]R": "P62805[aK@9][aK@13][aK@17](HIST1H4A)",
        "G[+56]K[+42]GGK[+42]GLGK[+42]GGAK[+56]R": "P62805[aK@6][aK@9][aK@13](HIST1H4A)",
        "G[+56]K[+42]GGK[+42]GLGK[+56]GGAK[+42]R": "P62805[aK@6][aK@9][aK@17](HIST1H4A)",
        "G[+56]K[+42]GGK[+42]GLGK[+42]GGAK[+42]R": "P62805[aK@6][aK@9][aK@13][aK@17](HIST1H4A)",
        "G[+56]K[+56]GGK[+56]GLGK[+56]GGAK[+70]R": "P62805[meK@17](HIST1H4A)",
        "K[+112.1]VLR": "P62805[me0K@21](HIST1H4A)",
        "K[+126.1]VLR": "P62805[meK@21](HIST1H4A)",
        "K[+84.1]VLR": "P62805[me2K@21](HIST1H4A)",
        "K[+98.1]VLR": "P62805[me3K@21](HIST1H4A)",
        "D[+56]AVTYTEHAK[+56]R": "P62805-NORM(HIST1H4A)",
        "Y[+56]QK[+28]STELLIR": "P68431[me2K@57](HIST1H3A)"
    }
    meta_data = pd.read_csv('static/data/P100-all-plates-metaData.csv', sep=',', header=None)
    #print meta_data
    for i in range(0, 96):
        # print str(meta_data[6][29 + i])
        # print petide2ptmAndGene[str(meta_data[6][29 + i])]
        # print str(meta_data[11][29 + i])
        # print str(meta_data[2][29 + i])
        row_peptide.append('Peptide: ' + str(meta_data[6][29 + i]) + "/" + peptide2ptmAndGene[str(meta_data[6][29 + i])])
        shorthand = peptide2ptmAndGene[str(meta_data[6][29 + i])]
        row_shorthand.append('PTM Proteins: ' + shorthand)
        row_protein.append('Protein: ' + str(meta_data[11][29 + i]))
        row_gene.append('Gene: ' + str(meta_data[2][29 + i]))
        print i

    col_cell_line = []
    col_pert_dose = []
    col_pert_time = []
    col_pert_name = []
    col_pert_signature = []

    col_pert_type = []
    first_loop = True




    # for i in range(0, 96):
    #     row_peptide.append('Peptide: ' + str(P100_data[6][29 + i]))
    #     shorthand = petide2ptmAndGene[str(P100_data[6][29 + i])]
    #     row_shorthand.append('PTM Proteins: ' + shorthand)
    #     row_protein.append('Protein: ' + str(P100_data[11][29 + i]))
    #     row_gene.append('Gene: ' + str(P100_data[2][29 + i]))

    # open('static/data/P100_processed_perturb_for_clustergramm_slicing.json') as f:

    with open('static/data/P100_processed_perturb_for_clustergramm_slicing.json') as f:
        data = json.load(f)
        cp_iter = 1
        for var in compounds:
            cp = str(var).upper()
            print cp
            if cp in data:
                for iter in data[cp]:

                    print cp_iter
                    col_cell_line.append('Cell Line: ' + iter["Cell_line"])
                    col_pert_dose.append('Dose: ' + iter["Dose"])
                    col_pert_time.append('Time: ' + iter["Time"])
                    col_pert_name.append(
                        iter["Perturbations"] + "/Dose/" + iter["Dose"] + "/Time/" + iter["Time"] + "/CellLine/" + iter[
                            "Cell_line"] + "/Signature ID/" + iter["Sig_id"])
                    col_pert_signature.append('Signature ID: ' + iter["Sig_id"])
                    col_pert_type.append('Type: ' + iter["Type"])

                    x1 = np.array(iter["data"])

                    # print "end of loop1"
                    x2 = x1[:, np.newaxis]
                    # print "end of loop2"
                    my_data_final = np.asarray(x2)
                    # print "end of loop3"
                    if first_loop:
                        my_data_con = my_data_final
                        first_loop = False

                    else:
                        my_data_con = np.hstack((my_data_con, my_data_final))


        # for i in range(0, 96):
        #     row_peptide.append('Peptide: ' + str(P100_data[6][29 + i]))
        #     shorthand = petide2ptmAndGene[str(P100_data[6][29 + i])]
        #     row_shorthand.append('PTM Proteins: ' + shorthand)
        #     row_protein.append('Protein: ' + str(P100_data[11][29 + i]))
        #     row_gene.append('Gene: ' + str(P100_data[2][29 + i]))

        arrays_rows2 = [np.array(row_peptide)]
        arrays_columns2 = [np.array(col_pert_name)]
        # print arrays_columns2
        # arrays_rows2 = np.load("static/data/P100_rowMetaDat")

        tuples_rows2 = list(zip(*arrays_rows2))
        # print tuples_rows2
        tuples_columns2 = list(zip(*arrays_columns2))
        # print tuples_columns2
        rows_labels2 = pd.MultiIndex.from_tuples(tuples_rows2)
        columns_labels2 = pd.MultiIndex.from_tuples(tuples_columns2)

        # my_data_con = my_data_con.fillna(0)
        print "after na to zero"
        a2 = np.nan_to_num(my_data_con)
        # a2 = my_data_con
        aa = a2.astype(np.float)



        print "a2.shape"
        print a2.shape
        print "columns_labels2.shape"
        print columns_labels2.shape
        print "rows_labels2.shape"
        print rows_labels2.shape
        # {'nop': row1, 'o0p': row2, 'zaz': row3, 'zax': row4, 'oof': row5, 'oye': row6}
        # df = pd.DataFrame(a, index=rows_labels, columns=columns_labels)
        df2 = pd.DataFrame(aa, index=rows_labels2, columns=columns_labels2)
        df2.to_csv("static/data/staurosporine.csv")
        jsonRet = df2.to_dict()

        print jsonRet
        return df2.to_json(orient='split')




@app.route("/api/clust/GCP/aggregatedforSlicing/<input_list>")
def GCPSlicing(input_list):
    compounds = []
    results = {}
    for var in input_list.split(","):
        compounds.append(var)

    #
    row_peptide = []
    row_protein = []
    row_gene = []
    row_shorthand = []
    peptide2ptmAndGene = {
        "IYQY[+80]IQSR": "Q13627[Y+80@321](DYRK1A)",
        "TPKDS[+80]PGIPPSANAHQLFR": "P51812[S+80@369](RPS6KA3)",
        "RNS[+80]SEASSGDFLDLK": "Q9UK76[S+80@87](JPT1)",
        "LPLVPES[+80]PRR": "Q86WB0[S+80@321](ZC3HC1)",
        "ANAS[+80]PQKPLDLK": "Q9Y618[S+80@956](NCOR2)",
        "LENS[+80]PLGEALR": "Q9NX40[S+80@108](OCIAD1)",
        "ANS[+80]FVGTAQYVSPELLTEK": "O15530[S+80@241](PDPK1)",
        "TNPPTQKPPS[+80]PPMSGR": "Q8IZP0[S+80@183](ABI1)",
        "SNS[+80]LPHSAVSNAGSK": "Q8TBZ3[S+80@434](WDR20)",
        "VGS[+80]LDNVGHLPAGGAVK": "P27816[S+80@1073](MAP4)",
        "AAPEAS[+80]SPPASPLQHLLPGK": "Q96TA1[S+80@691](FAM129B)",
        "S[+122]DKPDM[+16]AEIEKFDK": "P62328[S+122@2][M+16@7](TMSB4X)",
        "S[+122]DKPDMAEIEKFDK": "P62328[S+122@2](TMSB4X)",
        "SLS[+80]LGDKEISR": "Q9UMZ2[S+80@1075](SYNRG)",
        "DLVQPDKPAS[+80]PK": "Q6PJT7[S+80@515](ZC3H14)",
        "SPS[+80]PAHLPDDPKVAEK": "Q92615[S+80@601](LARP4B)",
        "S[+80]IQDLTVTGTEPGQVSSR": "O43318[S+80@439](MAP3K7)",
        "IHS[+80]PIIR": "O60885[S+80@1117](BRD4)",
        "TFS[+80]LTEVR": "O95239[S+80@801](KIF4A)",
        "SLVGS[+80]WLK": "Q6ICG6[S+80@362](KIAA0930)",
        "S[+80]PPAPGLQPMR": "P15408[S+80@200](FOSL2)",
        "LAS[+80]PELER": "P17535[S+80@100](JUND)",
        "IGPLGLS[+80]PK": "P30050[S+80@38](RPL12)",
        "TPS[+80]IQPSLLPHAAPFAK": "P35658[S+80@1023](NUP214)",
        "HAS[+80]PILPITEFSDIPR": "P42167[S+80@306](TMPO)",
        "LIPGPLS[+80]PVAR": "P48634[S+80@1219](PRRC2A)",
        "LGM[+16]LS[+80]PEGTC[+57]K": "P49327[S+80@207][M+16@205][C+57@212](FASN)",
        "LGMLS[+80]PEGTC[+57]K": "P49327[S+80@207][C+57@212](FASN)",
        "ISNLS[+80]PEEEQGLWK": "Q5HYJ3[S+80@193](FAM76B)",
        "VSMPDVELNLKS[+80]PK": "Q09666[S+80@3426](AHNAK)",
        "S[+122]DNGELEDKPPAPPVR": "Q13177[S+122@2](PAK2)",
        "KAYS[+80]FC[+57]GTVEYM[+16]APEVVNR": "Q15418[S+80@221][M+16@229][C+57@223](RPS6KA1)",
        "KAYS[+80]FC[+57]GTVEYMAPEVVNR": "Q15418[S+80@221][C+57@223](RPS6KA1)",
        "NDS[+80]WGSFDLR": "Q7Z417[S+80@652](NUFIP2)",
        "LEVTEIVKPS[+80]PK": "Q7Z6E9[S+80@1179](RBBP6)",
        "YGS[+80]PPQRDPNWNGER": "O15234[S+80@265](CASC3)",
        "QDDS[+80]PPRPIIGPALPPGFIK": "Q8IXQ4[S+80@105](GPALPP1)",
        "SFS[+80]ADNFIGIQR": "Q8N7R7[S+80@344](CCNYL1)",
        "VLS[+80]PLIIK": "Q8NCN4[S+80@403](RNF169)",
        "AGS[+80]PDVLR": "Q8NDX6[S+80@44](ZNF740)",
        "LGPGRPLPTFPTSEC[+57]TS[+80]DVEPDTR": "Q8TDD1[S+80@75][C+57@73](DDX54)",
        "LAAPSVSHVS[+80]PR": "Q8WXE1[S+80@224](ATRIP)",
        "VDDDS[+80]LGEFPVTNSR": "Q92785[S+80@142](DPF2)",
        "NEEPVRS[+80]PERR": "Q92922[S+80@310](SMARCC1)",
        "LFIIRGS[+80]PQQIDHAK": "Q92945[S+80@480](KHSRP)",
        "S[+80]IEVENDFLPVEK": "Q96B97[S+80@230](SH3KBP1)",
        "TAPTLS[+80]PEHWK": "Q96JM3[S+80@405](CHAMP1)",
        "VLS[+80]PTAAKPSPFEGK": "Q96QC0[S+80@313](PPP1R10)",
        "SSDQPLTVPVS[+80]PK": "Q9ULW0[S+80@738](TPX2)",
        "FYETKEESYS[+80]PSKDR": "Q96T23[S+80@473](RSF1)",
        "SDS[+80]PENKYSDSTGHSK": "Q9BTA9[S+80@64](WAC)",
        "S[+80]IPLSIK": "Q9C0C9[S+80@515](UBE2O)",
        "RLS[+80]QSDEDVIR": "Q9H7D7[S+80@121](WDR26)",
        "ATS[+80]PVKSTTSITDAK": "Q9NQW6[S+80@295](ANLN)",
        "ALGS[+80]PTKQLLPC[+57]EMAC[+57]NEK": "Q9NR45[S+80@275][C+57@283][C+57@287](NANS)",
        "YLLGDAPVS[+80]PSSQK": "Q9NYB0[S+80@203](TERF2IP)",
        "ANS[+80]PEKPPEAGAAHKPR": "Q9UFC0[S+80@212](LRWD1)",
        "SEVQQPVHPKPLS[+80]PDSR": "Q9UHB6[S+80@362](LIMA1)",
        "ETPHS[+80]PGVEDAPIAK": "Q9UHB6[S+80@490](LIMA1)",
        "SQS[+80]PHYFR": "Q9UKJ3[S+80@1035](GPATCH8)",
        "DRS[+80]SPPPGYIPDELHQVAR": "Q9Y2U5[S+80@163](MAP3K2)",
        "SPALKS[+80]PLQSVVVR": "Q9Y2W1[S+80@253](THRAP3)",
        "AFGSGIDIKPGT[+80]PPIAGR": "Q9Y520[T+80@2673](PRRC2C)",
        "SFS[+80]SQRPVDR": "Q9Y520[S+80@1544](PRRC2C)",
        "VYT[+80]HEVVTLWYR": "P06493[T+80@161](CDK1)",
        "SST[+80]PLPTISSSAENTR": "P42167[T+80@160](TMPO)",
        "QIT[+80]MEELVR": "Q15149[T+80@4030](PLEC)",
        "TQLWASEPGT[+80]PPLPTSLPSQNPILK": "Q9BXP5[T+80@544](SRRT)",
        "ALPQT[+80]PRPR": "Q9UQ35[T+80@1492](SRRM2)",
        "SMS[+80]VDLSHIPLKDPLLFK": "A0JNW5[S+80@935](UHRF1BP1L)",
        "S[+80]PTGPSNSFLANMGGTVAHK": "Q96I25[S+80@222](RBM17)",
        "S[+80]LTAHSLLPLAEK": "Q86VI3[S+80@1424](IQGAP3)",
        "S[+80]FAGNLNTYKR": "Q01813[S+80@386](PFKP)",
        "HRPS[+80]PPATPPPK": "Q8IYB3[S+80@402](SRRM1)",
        "LHS[+80]APNLSDLHVVRPK": "O75385[S+80@556](ULK1)",
        "TLGRRDS[+80]SDDWEIPDGQITVGQR": "P15056[S+80@446](BRAF)",
        "A[+42]TTATM[+16]ATSGS[+80]AR": "P38919[S+80@12][M+16@7][A+42@2](EIF4A3)",
        "A[+42]TTATMATSGS[+80]AR": "P38919[S+80@12][A+42@2](EIF4A3)",
        "IHVSRS[+80]PTRPR": "Q499Z4[S+80@189](ZNF672)",
        "RPHS[+80]PEKAFSSNPVVR": "Q53F19[S+80@500](NCBP3)",
        "KPNIFYSGPAS[+80]PARPR": "Q6PL18[S+80@327](ATAD2)",
        "TEFLDLDNSPLSPPS[+80]PR": "Q8NCF5[S+80@204](NFATC2IP)",
        "QGSGRES[+80]PSLASR": "Q8WWM7[S+80@339](ATXN2L)",
        # "TQLWASEPGT[+80]PPLPTSLPSQNPILK": "Q8WWM7[S+80@339](SRRM2)",
        "LQS[+80]EPESIR": "P09496[S+80@105](CLTA)",
        "RLIS[+80]PYKK": "O14929[S+80@361](HAT1)",
        "LLEDS[+80]EESSEETVSR": "O60231[S+80@103](DHX16)",
        "S[+80]PPAPGLQPM[+16]R": "P15408[S+80@200][M+16@209](FOSL2)",
        "RRLS[+80]SLR": "P62753[S+80@235](RPS6)",
        "RLS[+80]ESQLSFRR": "Q96PK6[S+80@618](RBM14)",
        "RLS[+80]LPGLLSQVSPR": "Q96Q42[S+80@483](ALS2)",
        "SPDKPGGS[+80]PSASRR": "Q9Y3T9[S+80@56](NOC2L)",
        "HLPS[+80]PPTLDSIITEYLR": "Q9Y4B6[S+80@1000](DCAF1)",
        "ST[+80]FHAGQLR": "Q7KZI7[T+80@596](MARK2)",
        "S[+80]LTNSHLEKK": "Q9H2H9[S+80@52](SLC38A1)",
        "LQTPNT[+80]FPKR": "Q14978[T+80@610](NOLC1)",
        "QIT[+80]M[+16]EELVR": "Q15149[M+16@4031][T+80@4030](PLEC)",

        "T[+56]K[+56]QTAR": "P68431[me0K@5](HIST1H3A)",
        "T[+56]K[+70]QTAR": "P68431[meK@5](HIST1H3A)",
        "T[+56]K[+28]QTAR": "P68431[me2K@5](HIST1H3A)",
        "T[+56]K[+42]QTAR-me3K": "P68431[me3K@5](HIST1H3A)",
        "T[+56]K[+42]QTAR-aK": "P68431[aK@5](HIST1H3A)",
        "K[+112.1]STGGK[+56]APR": "P68431[me0K@10][a0K@15](HIST1H3A)",
        "K[+126.1]STGGK[+56]APR": "P68431[meK@10](HIST1H3A)",
        "K[+84.1]STGGK[+56]APR": "P68431[me2K@10](HIST1H3A)",
        "K[+98.1]STGGK[+56]APR": "P68431[me3K@10](HIST1H3A)",
        "K[+98]STGGK[+56]APR": "P68431[aK@10](HIST1H3A)",
        "K[+112.1]STGGK[+42]APR": "P68431[aK@15](HIST1H3A)",
        "K[+126.1]STGGK[+42]APR": "P68431[meK@10][aK@15](HIST1H3A)",
        "K[+84.1]STGGK[+42]APR": "P68431[me2K@10][aK@15](HIST1H3A)",
        "K[+98.1]STGGK[+42]APR": "P68431[me3K@10][aK@15](HIST1H3A)",
        "K[+98]STGGK[+42]APR": "P68431[aK@10][aK@15](HIST1H3A)",
        "K[+112.1]S[+80]TGGK[+56]APR": "P68431[pS@11](HIST1H3A)",
        "K[+126.1]S[+80]TGGK[+56]APR": "P68431[meK@10][pS@11](HIST1H3A)",
        "K[+84.1]S[+80]TGGK[+56]APR": "P68431[me2K@10][pS@11](HIST1H3A)",
        "K[+98.1]S[+80]TGGK[+56]APR": "P68431[me3K@10][pS@11](HIST1H3A)",
        "K[+98]S[+80]TGGK[+56]APR": "P68431[aK@10][pS@11](HIST1H3A)",
        "K[+112.1]S[+80]TGGK[+42]APR": "P68431[pS@11][aK@15](HIST1H3A)",
        "K[+126.1]S[+80]TGGK[+42]APR": "P68431[meK@10][pS@11][aK@15](HIST1H3A)",
        "K[+84.1]S[+80]TGGK[+42]APR": "P68431[me2K@10][pS@11][aK@15](HIST1H3A)",
        "K[+98.1]S[+80]TGGK[+42]APR": "P68431[me3K@10][pS@11][aK@15](HIST1H3A)",
        "K[+98]S[+80]TGGK[+42]APR": "P68431[aK@10][pS@11][aK@15](HIST1H3A)",

        "K[+112.1]QLATK[+56]AAR": "P68431[a0K@19][a0K@24](HIST1H3A)",
        "K[+98]QLATK[+56]AAR": "P68431[aK@19](HIST1H3A)",
        "K[+112.1]QLATK[+42]AAR": "P68431[aK@24](HIST1H3A)",
        "K[+98]QLATK[+42]AAR": "P68431[aK@19][aK@24](HIST1H3A)",

        "K[+226.1]QLATK[+56]AAR": "P68431[ubK@19][a0K@24](HIST1H3A)",
        "K[+112.1]QLATK[+170.1]AAR": "P68431[ubK@24](HIST1H3A)",
        "K[+112.1]SAPATGGVK[+56]K[+56]PHR": "P68431[me0K@28][me0K@37](HIST1H3A)",
        "K[+126.1]SAPATGGVK[+56]K[+56]PHR": "P68431[meK@28][me0K@37](HIST1H3A)",
        "K[+126.1]SAPATGGVK[+70]K[+56]PHR": "P68431[meK@28][meK@37](HIST1H3A)",
        "K[+126.1]SAPATGGVK[+28]K[+56]PHR": "P68431[meK@28][me2K@37](HIST1H3A)",
        "K[+126.1]SAPATGGVK[+42]K[+56]PHR": "P68431[meK@28][me3K@37](HIST1H3A)",
        "K[+84.1]SAPATGGVK[+56]K[+56]PHR": "P68431[me2K@28][me0K@37](HIST1H3A)",
        "K[+84.1]SAPATGGVK[+70]K[+56]PHR": "P68431[me2K@28][meK@37](HIST1H3A)",
        "K[+84.1]SAPATGGVK[+28]K[+56]PHR": "P68431[me2K@28][me2K@37](HIST1H3A)",
        "K[+84.1]SAPATGGVK[+42]K[+56]PHR": "P68431[me2K@28][me3K@37](HIST1H3A)",
        "K[+98.1]SAPATGGVK[+56]K[+56]PHR": "P68431[me3K@28][me0K@37](HIST1H3A)",
        "K[+98.1]SAPATGGVK[+70]K[+56]PHR": "P68431[me3K@28][meK@37](HIST1H3A)",
        "K[+98.1]SAPATGGVK[+28]K[+56]PHR": "P68431[me3K@28][me2K@37](HIST1H3A)",
        "K[+98.1]SAPATGGVK[+42]K[+56]PHR": "P68431[me3K@28][me3K@37](HIST1H3A)",
        "K[+98]SAPATGGVK[+56]K[+56]PHR": "P68431[aK@28][me0K@37](HIST1H3A)",
        "K[+98]SAPATGGVK[+70]K[+56]PHR": "P68431[aK@28][meK@37](HIST1H3A)",
        "K[+98]SAPATGGVK[+28]K[+56]PHR": "P68431[aK@28][me2K@37](HIST1H3A)",
        "K[+98]SAPATGGVK[+42]K[+56]PHR": "P68431[aK@28][me3K@37](HIST1H3A)",
        "K[+112.1]SAPSTGGVK[+56]K[+56]PHR": "P84243[me0K@28][me0K@37](H3F3A)",
        "Y[+56]RPGTVALR": "P68431-NORM(HIST1H3A)",
        "Y[+56]QK[+56]STELLIR": "P68431[me0K@57](HIST1H3A)",
        "E[+56]IAQDFK[+56]TDLR": "P68431[me0K@80](HIST1H3A)",
        "E[+56]IAQDFK[+70]TDLR": "P68431[meK@80](HIST1H3A)",
        "E[+56]IAQDFK[+28]TDLR": "P68431[me2K@80](HIST1H3A)",
        "K[+112.1]SAPATGGVK[+70]K[+56]PHR": "P68431[me0K@28][meK@37](HIST1H3A)",
        "K[+112.1]SAPATGGVK[+28]K[+56]PHR": "P68431[me0K@28][me2K@37](HIST1H3A)",
        "K[+112.1]SAPATGGVK[+42]K[+56]PHR": "P68431[me0K@28][me3K@37](HIST1H3A)",
        "Y[+56]QK[+42]STELLIR": "P68431[aK@57](HIST1H3A)",
        "G[+56]K[+56]GGK[+56]GLGK[+56]GGAK[+56]R": "P62805[a0K@6][a0K@13][a0K@17](HIST1H4A)",
        "G[+56]K[+42]GGK[+56]GLGK[+56]GGAK[+56]R": "P62805[aK@6](HIST1H4A)",
        "G[+56]K[+56]GGK[+56]GLGK[+42]GGAK[+56]R": "P62805[aK@13](HIST1H4A)",
        "G[+56]K[+56]GGK[+56]GLGK[+56]GGAK[+42]R": "P62805[aK@17](HIST1H4A)",
        "G[+56]K[+56]GGK[+42]GLGK[+42]GGAK[+56]R": "P62805[aK@9][aK@13](HIST1H4A)",
        "G[+56]K[+42]GGK[+42]GLGK[+56]GGAK[+56]R": "P62805[aK@6][aK@9](HIST1H4A)",
        "G[+56]K[+42]GGK[+56]GLGK[+56]GGAK[+42]R": "P62805[aK@6][aK@17](HIST1H4A)",
        "G[+56]K[+56]GGK[+56]GLGK[+42]GGAK[+42]R": "P62805[aK@13][aK@17](HIST1H4A)",
        "G[+56]K[+56]GGK[+42]GLGK[+42]GGAK[+42]R": "P62805[aK@9][aK@13][aK@17](HIST1H4A)",
        "G[+56]K[+42]GGK[+42]GLGK[+42]GGAK[+56]R": "P62805[aK@6][aK@9][aK@13](HIST1H4A)",
        "G[+56]K[+42]GGK[+42]GLGK[+56]GGAK[+42]R": "P62805[aK@6][aK@9][aK@17](HIST1H4A)",
        "G[+56]K[+42]GGK[+42]GLGK[+42]GGAK[+42]R": "P62805[aK@6][aK@9][aK@13][aK@17](HIST1H4A)",
        "G[+56]K[+56]GGK[+56]GLGK[+56]GGAK[+70]R": "P62805[meK@17](HIST1H4A)",
        "K[+112.1]VLR": "P62805[me0K@21](HIST1H4A)",
        "K[+126.1]VLR": "P62805[meK@21](HIST1H4A)",
        "K[+84.1]VLR": "P62805[me2K@21](HIST1H4A)",
        "K[+98.1]VLR": "P62805[me3K@21](HIST1H4A)",
        "D[+56]AVTYTEHAK[+56]R": "P62805-NORM(HIST1H4A)",
        "Y[+56]QK[+28]STELLIR": "P68431[me2K@57](HIST1H3A)"
    }
    meta_data = pd.read_csv('static/data/GCP-all-plates-metaData.csv', sep=',', header=None)
    #print meta_data

    for i in range(0, 79):

        peptide = str(meta_data[4][28 + i])
        if(peptide == 'T[+56]K[+42]QTAR'):
            print str(meta_data[2][28 + i])
            if(str(meta_data[2][28 + i]) == 'BI10006'):

                peptide = "T[+56]K[+42]QTAR-me3K"
            else:
                peptide = "T[+56]K[+42]QTAR-aK"
        print peptide
        row_peptide.append('Peptide: ' + peptide)

        shorthand = peptide2ptmAndGene[peptide]

        row_shorthand.append('PTM Proteins: ' + shorthand)
        row_protein.append('Protein: ' + str(meta_data[9][28 + i]))
        row_gene.append('Gene: ' + str(meta_data[6][28 + i]))



    col_cell_line = []
    col_pert_dose = []
    col_pert_time = []
    col_pert_name = []
    col_pert_signature = []

    col_pert_type = []
    first_loop = True




    # for i in range(0, 96):
    #     row_peptide.append('Peptide: ' + str(P100_data[6][29 + i]))
    #     shorthand = petide2ptmAndGene[str(P100_data[6][29 + i])]
    #     row_shorthand.append('PTM Proteins: ' + shorthand)
    #     row_protein.append('Protein: ' + str(P100_data[11][29 + i]))
    #     row_gene.append('Gene: ' + str(P100_data[2][29 + i]))

    # open('static/data/P100_processed_perturb_for_clustergramm_slicing.json') as f:



    with open('static/data/GCP_processed_perturb_for_clustergramm_slicing.json') as f:
        data = json.load(f)
        cp_iter = 1
        for var in compounds:
            cp = str(var).upper()
            print cp

            if cp in data:
                for iter in data[cp]:

                    print cp_iter
                    cp_iter += 1
                    col_cell_line.append('Cell Line: ' + iter["Cell_line"])
                    col_pert_dose.append('Dose: ' + iter["Dose"])
                    col_pert_time.append('Time: ' + iter["Time"])
                    col_pert_name.append('Perturbations: ' + iter["Perturbations"])
                    col_pert_signature.append('Signature ID: ' + iter["Sig_id"])
                    col_pert_type.append('Type: ' + iter["Type"])

                    x1 = np.array(iter["data"])

                    # print "end of loop1"
                    x2 = x1[:, np.newaxis]
                    # print "end of loop2"
                    my_data_final = np.asarray(x2)
                    # print "end of loop3"
                    if first_loop:
                        my_data_con = my_data_final
                        first_loop = False

                    else:
                        my_data_con = np.hstack((my_data_con, my_data_final))


        # for i in range(0, 96):
        #     row_peptide.append('Peptide: ' + str(P100_data[6][29 + i]))
        #     shorthand = petide2ptmAndGene[str(P100_data[6][29 + i])]
        #     row_shorthand.append('PTM Proteins: ' + shorthand)
        #     row_protein.append('Protein: ' + str(P100_data[11][29 + i]))
        #     row_gene.append('Gene: ' + str(P100_data[2][29 + i]))

        arrays_rows2 = [np.array(row_shorthand), np.array(row_peptide), np.array(row_protein), np.array(row_gene)]
        arrays_columns2 = [np.array(col_pert_name), np.array(col_pert_signature), np.array(col_pert_time),
                           np.array(col_pert_dose), np.array(col_pert_type), np.array(col_cell_line)]
        # print arrays_columns2
        # arrays_rows2 = np.load("static/data/P100_rowMetaDat")

        tuples_rows2 = list(zip(*arrays_rows2))
        # print tuples_rows2
        tuples_columns2 = list(zip(*arrays_columns2))
        # print tuples_columns2
        rows_labels2 = pd.MultiIndex.from_tuples(tuples_rows2)
        columns_labels2 = pd.MultiIndex.from_tuples(tuples_columns2)

        # my_data_con = my_data_con.fillna(0)
        print "after na to zero"
        a2 = np.nan_to_num(my_data_con)
        # a2 = my_data_con
        aa = a2.astype(np.float)

        print "a2.shape"
        print a2.shape
        print "columns_labels2.shape"
        print columns_labels2.shape
        print "rows_labels2.shape"
        print rows_labels2.shape
        # {'nop': row1, 'o0p': row2, 'zaz': row3, 'zax': row4, 'oof': row5, 'oye': row6}
        # df = pd.DataFrame(a, index=rows_labels, columns=columns_labels)
        df2 = pd.DataFrame(aa, index=rows_labels2, columns=columns_labels2)



            #print results
    print "finished"
    GCP_all_json = loadFile.make_json_from_txt(df2)
    return GCP_all_json





@app.route("/api/clust/P100/aggregatedforSlicing/<input_list>")
def P100Slicing(input_list):
    compounds = []
    results = {}
    for var in input_list.split(","):
        compounds.append(var)

    #
    row_peptide = []
    row_protein = []
    row_gene = []
    row_shorthand = []
    peptide2ptmAndGene = {
        "IYQY[+80]IQSR": "Q13627[Y+80@321](DYRK1A)",
        "TPKDS[+80]PGIPPSANAHQLFR": "P51812[S+80@369](RPS6KA3)",
        "RNS[+80]SEASSGDFLDLK": "Q9UK76[S+80@87](JPT1)",
        "LPLVPES[+80]PRR": "Q86WB0[S+80@321](ZC3HC1)",
        "ANAS[+80]PQKPLDLK": "Q9Y618[S+80@956](NCOR2)",
        "LENS[+80]PLGEALR": "Q9NX40[S+80@108](OCIAD1)",
        "ANS[+80]FVGTAQYVSPELLTEK": "O15530[S+80@241](PDPK1)",
        "TNPPTQKPPS[+80]PPMSGR": "Q8IZP0[S+80@183](ABI1)",
        "SNS[+80]LPHSAVSNAGSK": "Q8TBZ3[S+80@434](WDR20)",
        "VGS[+80]LDNVGHLPAGGAVK": "P27816[S+80@1073](MAP4)",
        "AAPEAS[+80]SPPASPLQHLLPGK": "Q96TA1[S+80@691](FAM129B)",
        "S[+122]DKPDM[+16]AEIEKFDK": "P62328[S+122@2][M+16@7](TMSB4X)",
        "S[+122]DKPDMAEIEKFDK": "P62328[S+122@2](TMSB4X)",
        "SLS[+80]LGDKEISR": "Q9UMZ2[S+80@1075](SYNRG)",
        "DLVQPDKPAS[+80]PK": "Q6PJT7[S+80@515](ZC3H14)",
        "SPS[+80]PAHLPDDPKVAEK": "Q92615[S+80@601](LARP4B)",
        "S[+80]IQDLTVTGTEPGQVSSR": "O43318[S+80@439](MAP3K7)",
        "IHS[+80]PIIR": "O60885[S+80@1117](BRD4)",
        "TFS[+80]LTEVR": "O95239[S+80@801](KIF4A)",
        "SLVGS[+80]WLK": "Q6ICG6[S+80@362](KIAA0930)",
        "S[+80]PPAPGLQPMR": "P15408[S+80@200](FOSL2)",
        "LAS[+80]PELER": "P17535[S+80@100](JUND)",
        "IGPLGLS[+80]PK": "P30050[S+80@38](RPL12)",
        "TPS[+80]IQPSLLPHAAPFAK": "P35658[S+80@1023](NUP214)",
        "HAS[+80]PILPITEFSDIPR": "P42167[S+80@306](TMPO)",
        "LIPGPLS[+80]PVAR": "P48634[S+80@1219](PRRC2A)",
        "LGM[+16]LS[+80]PEGTC[+57]K": "P49327[S+80@207][M+16@205][C+57@212](FASN)",
        "LGMLS[+80]PEGTC[+57]K": "P49327[S+80@207][C+57@212](FASN)",
        "ISNLS[+80]PEEEQGLWK": "Q5HYJ3[S+80@193](FAM76B)",
        "VSMPDVELNLKS[+80]PK": "Q09666[S+80@3426](AHNAK)",
        "S[+122]DNGELEDKPPAPPVR": "Q13177[S+122@2](PAK2)",
        "KAYS[+80]FC[+57]GTVEYM[+16]APEVVNR": "Q15418[S+80@221][M+16@229][C+57@223](RPS6KA1)",
        "KAYS[+80]FC[+57]GTVEYMAPEVVNR": "Q15418[S+80@221][C+57@223](RPS6KA1)",
        "NDS[+80]WGSFDLR": "Q7Z417[S+80@652](NUFIP2)",
        "LEVTEIVKPS[+80]PK": "Q7Z6E9[S+80@1179](RBBP6)",
        "YGS[+80]PPQRDPNWNGER": "O15234[S+80@265](CASC3)",
        "QDDS[+80]PPRPIIGPALPPGFIK": "Q8IXQ4[S+80@105](GPALPP1)",
        "SFS[+80]ADNFIGIQR": "Q8N7R7[S+80@344](CCNYL1)",
        "VLS[+80]PLIIK": "Q8NCN4[S+80@403](RNF169)",
        "AGS[+80]PDVLR": "Q8NDX6[S+80@44](ZNF740)",
        "LGPGRPLPTFPTSEC[+57]TS[+80]DVEPDTR": "Q8TDD1[S+80@75][C+57@73](DDX54)",
        "LAAPSVSHVS[+80]PR": "Q8WXE1[S+80@224](ATRIP)",
        "VDDDS[+80]LGEFPVTNSR": "Q92785[S+80@142](DPF2)",
        "NEEPVRS[+80]PERR": "Q92922[S+80@310](SMARCC1)",
        "LFIIRGS[+80]PQQIDHAK": "Q92945[S+80@480](KHSRP)",
        "S[+80]IEVENDFLPVEK": "Q96B97[S+80@230](SH3KBP1)",
        "TAPTLS[+80]PEHWK": "Q96JM3[S+80@405](CHAMP1)",
        "VLS[+80]PTAAKPSPFEGK": "Q96QC0[S+80@313](PPP1R10)",
        "SSDQPLTVPVS[+80]PK": "Q9ULW0[S+80@738](TPX2)",
        "FYETKEESYS[+80]PSKDR": "Q96T23[S+80@473](RSF1)",
        "SDS[+80]PENKYSDSTGHSK": "Q9BTA9[S+80@64](WAC)",
        "S[+80]IPLSIK": "Q9C0C9[S+80@515](UBE2O)",
        "RLS[+80]QSDEDVIR": "Q9H7D7[S+80@121](WDR26)",
        "ATS[+80]PVKSTTSITDAK": "Q9NQW6[S+80@295](ANLN)",
        "ALGS[+80]PTKQLLPC[+57]EMAC[+57]NEK": "Q9NR45[S+80@275][C+57@283][C+57@287](NANS)",
        "YLLGDAPVS[+80]PSSQK": "Q9NYB0[S+80@203](TERF2IP)",
        "ANS[+80]PEKPPEAGAAHKPR": "Q9UFC0[S+80@212](LRWD1)",
        "SEVQQPVHPKPLS[+80]PDSR": "Q9UHB6[S+80@362](LIMA1)",
        "ETPHS[+80]PGVEDAPIAK": "Q9UHB6[S+80@490](LIMA1)",
        "SQS[+80]PHYFR": "Q9UKJ3[S+80@1035](GPATCH8)",
        "DRS[+80]SPPPGYIPDELHQVAR": "Q9Y2U5[S+80@163](MAP3K2)",
        "SPALKS[+80]PLQSVVVR": "Q9Y2W1[S+80@253](THRAP3)",
        "AFGSGIDIKPGT[+80]PPIAGR": "Q9Y520[T+80@2673](PRRC2C)",
        "SFS[+80]SQRPVDR": "Q9Y520[S+80@1544](PRRC2C)",
        "VYT[+80]HEVVTLWYR": "P06493[T+80@161](CDK1)",
        "SST[+80]PLPTISSSAENTR": "P42167[T+80@160](TMPO)",
        "QIT[+80]MEELVR": "Q15149[T+80@4030](PLEC)",
        "TQLWASEPGT[+80]PPLPTSLPSQNPILK": "Q9BXP5[T+80@544](SRRT)",
        "ALPQT[+80]PRPR": "Q9UQ35[T+80@1492](SRRM2)",
        "SMS[+80]VDLSHIPLKDPLLFK": "A0JNW5[S+80@935](UHRF1BP1L)",
        "S[+80]PTGPSNSFLANMGGTVAHK": "Q96I25[S+80@222](RBM17)",
        "S[+80]LTAHSLLPLAEK": "Q86VI3[S+80@1424](IQGAP3)",
        "S[+80]FAGNLNTYKR": "Q01813[S+80@386](PFKP)",
        "HRPS[+80]PPATPPPK": "Q8IYB3[S+80@402](SRRM1)",
        "LHS[+80]APNLSDLHVVRPK": "O75385[S+80@556](ULK1)",
        "TLGRRDS[+80]SDDWEIPDGQITVGQR": "P15056[S+80@446](BRAF)",
        "A[+42]TTATM[+16]ATSGS[+80]AR": "P38919[S+80@12][M+16@7][A+42@2](EIF4A3)",
        "A[+42]TTATMATSGS[+80]AR": "P38919[S+80@12][A+42@2](EIF4A3)",
        "IHVSRS[+80]PTRPR": "Q499Z4[S+80@189](ZNF672)",
        "RPHS[+80]PEKAFSSNPVVR": "Q53F19[S+80@500](NCBP3)",
        "KPNIFYSGPAS[+80]PARPR": "Q6PL18[S+80@327](ATAD2)",
        "TEFLDLDNSPLSPPS[+80]PR": "Q8NCF5[S+80@204](NFATC2IP)",
        "QGSGRES[+80]PSLASR": "Q8WWM7[S+80@339](ATXN2L)",
        # "TQLWASEPGT[+80]PPLPTSLPSQNPILK": "Q8WWM7[S+80@339](SRRM2)",
        "LQS[+80]EPESIR": "P09496[S+80@105](CLTA)",
        "RLIS[+80]PYKK": "O14929[S+80@361](HAT1)",
        "LLEDS[+80]EESSEETVSR": "O60231[S+80@103](DHX16)",
        "S[+80]PPAPGLQPM[+16]R": "P15408[S+80@200][M+16@209](FOSL2)",
        "RRLS[+80]SLR": "P62753[S+80@235](RPS6)",
        "RLS[+80]ESQLSFRR": "Q96PK6[S+80@618](RBM14)",
        "RLS[+80]LPGLLSQVSPR": "Q96Q42[S+80@483](ALS2)",
        "SPDKPGGS[+80]PSASRR": "Q9Y3T9[S+80@56](NOC2L)",
        "HLPS[+80]PPTLDSIITEYLR": "Q9Y4B6[S+80@1000](DCAF1)",
        "ST[+80]FHAGQLR": "Q7KZI7[T+80@596](MARK2)",
        "S[+80]LTNSHLEKK": "Q9H2H9[S+80@52](SLC38A1)",
        "LQTPNT[+80]FPKR": "Q14978[T+80@610](NOLC1)",
        "QIT[+80]M[+16]EELVR": "Q15149[M+16@4031][T+80@4030](PLEC)",

        "T[+56]K[+56]QTAR": "P68431[me0K@5](HIST1H3A)",
        "T[+56]K[+70]QTAR": "P68431[meK@5](HIST1H3A)",
        "T[+56]K[+28]QTAR": "P68431[me2K@5](HIST1H3A)",
        "T[+56]K[+42]QTAR-me3K": "P68431[me3K@5](HIST1H3A)",
        "T[+56]K[+42]QTAR-aK": "P68431[aK@5](HIST1H3A)",
        "K[+112.1]STGGK[+56]APR": "P68431[me0K@10][a0K@15](HIST1H3A)",
        "K[+126.1]STGGK[+56]APR": "P68431[meK@10](HIST1H3A)",
        "K[+84.1]STGGK[+56]APR": "P68431[me2K@10](HIST1H3A)",
        "K[+98.1]STGGK[+56]APR": "P68431[me3K@10](HIST1H3A)",
        "K[+98]STGGK[+56]APR": "P68431[aK@10](HIST1H3A)",
        "K[+112.1]STGGK[+42]APR": "P68431[aK@15](HIST1H3A)",
        "K[+126.1]STGGK[+42]APR": "P68431[meK@10][aK@15](HIST1H3A)",
        "K[+84.1]STGGK[+42]APR": "P68431[me2K@10][aK@15](HIST1H3A)",
        "K[+98.1]STGGK[+42]APR": "P68431[me3K@10][aK@15](HIST1H3A)",
        "K[+98]STGGK[+42]APR": "P68431[aK@10][aK@15](HIST1H3A)",
        "K[+112.1]S[+80]TGGK[+56]APR": "P68431[pS@11](HIST1H3A)",
        "K[+126.1]S[+80]TGGK[+56]APR": "P68431[meK@10][pS@11](HIST1H3A)",
        "K[+84.1]S[+80]TGGK[+56]APR": "P68431[me2K@10][pS@11](HIST1H3A)",
        "K[+98.1]S[+80]TGGK[+56]APR": "P68431[me3K@10][pS@11](HIST1H3A)",
        "K[+98]S[+80]TGGK[+56]APR": "P68431[aK@10][pS@11](HIST1H3A)",
        "K[+112.1]S[+80]TGGK[+42]APR": "P68431[pS@11][aK@15](HIST1H3A)",
        "K[+126.1]S[+80]TGGK[+42]APR": "P68431[meK@10][pS@11][aK@15](HIST1H3A)",
        "K[+84.1]S[+80]TGGK[+42]APR": "P68431[me2K@10][pS@11][aK@15](HIST1H3A)",
        "K[+98.1]S[+80]TGGK[+42]APR": "P68431[me3K@10][pS@11][aK@15](HIST1H3A)",
        "K[+98]S[+80]TGGK[+42]APR": "P68431[aK@10][pS@11][aK@15](HIST1H3A)",

        "K[+112.1]QLATK[+56]AAR": "P68431[a0K@19][a0K@24](HIST1H3A)",
        "K[+98]QLATK[+56]AAR": "P68431[aK@19](HIST1H3A)",
        "K[+112.1]QLATK[+42]AAR": "P68431[aK@24](HIST1H3A)",
        "K[+98]QLATK[+42]AAR": "P68431[aK@19][aK@24](HIST1H3A)",

        "K[+226.1]QLATK[+56]AAR": "P68431[ubK@19][a0K@24](HIST1H3A)",
        "K[+112.1]QLATK[+170.1]AAR": "P68431[ubK@24](HIST1H3A)",
        "K[+112.1]SAPATGGVK[+56]K[+56]PHR": "P68431[me0K@28][me0K@37](HIST1H3A)",
        "K[+126.1]SAPATGGVK[+56]K[+56]PHR": "P68431[meK@28][me0K@37](HIST1H3A)",
        "K[+126.1]SAPATGGVK[+70]K[+56]PHR": "P68431[meK@28][meK@37](HIST1H3A)",
        "K[+126.1]SAPATGGVK[+28]K[+56]PHR": "P68431[meK@28][me2K@37](HIST1H3A)",
        "K[+126.1]SAPATGGVK[+42]K[+56]PHR": "P68431[meK@28][me3K@37](HIST1H3A)",
        "K[+84.1]SAPATGGVK[+56]K[+56]PHR": "P68431[me2K@28][me0K@37](HIST1H3A)",
        "K[+84.1]SAPATGGVK[+70]K[+56]PHR": "P68431[me2K@28][meK@37](HIST1H3A)",
        "K[+84.1]SAPATGGVK[+28]K[+56]PHR": "P68431[me2K@28][me2K@37](HIST1H3A)",
        "K[+84.1]SAPATGGVK[+42]K[+56]PHR": "P68431[me2K@28][me3K@37](HIST1H3A)",
        "K[+98.1]SAPATGGVK[+56]K[+56]PHR": "P68431[me3K@28][me0K@37](HIST1H3A)",
        "K[+98.1]SAPATGGVK[+70]K[+56]PHR": "P68431[me3K@28][meK@37](HIST1H3A)",
        "K[+98.1]SAPATGGVK[+28]K[+56]PHR": "P68431[me3K@28][me2K@37](HIST1H3A)",
        "K[+98.1]SAPATGGVK[+42]K[+56]PHR": "P68431[me3K@28][me3K@37](HIST1H3A)",
        "K[+98]SAPATGGVK[+56]K[+56]PHR": "P68431[aK@28][me0K@37](HIST1H3A)",
        "K[+98]SAPATGGVK[+70]K[+56]PHR": "P68431[aK@28][meK@37](HIST1H3A)",
        "K[+98]SAPATGGVK[+28]K[+56]PHR": "P68431[aK@28][me2K@37](HIST1H3A)",
        "K[+98]SAPATGGVK[+42]K[+56]PHR": "P68431[aK@28][me3K@37](HIST1H3A)",
        "K[+112.1]SAPSTGGVK[+56]K[+56]PHR": "P84243[me0K@28][me0K@37](H3F3A)",
        "Y[+56]RPGTVALR": "P68431-NORM(HIST1H3A)",
        "Y[+56]QK[+56]STELLIR": "P68431[me0K@57](HIST1H3A)",
        "E[+56]IAQDFK[+56]TDLR": "P68431[me0K@80](HIST1H3A)",
        "E[+56]IAQDFK[+70]TDLR": "P68431[meK@80](HIST1H3A)",
        "E[+56]IAQDFK[+28]TDLR": "P68431[me2K@80](HIST1H3A)",
        "K[+112.1]SAPATGGVK[+70]K[+56]PHR": "P68431[me0K@28][meK@37](HIST1H3A)",
        "K[+112.1]SAPATGGVK[+28]K[+56]PHR": "P68431[me0K@28][me2K@37](HIST1H3A)",
        "K[+112.1]SAPATGGVK[+42]K[+56]PHR": "P68431[me0K@28][me3K@37](HIST1H3A)",
        "Y[+56]QK[+42]STELLIR": "P68431[aK@57](HIST1H3A)",
        "G[+56]K[+56]GGK[+56]GLGK[+56]GGAK[+56]R": "P62805[a0K@6][a0K@13][a0K@17](HIST1H4A)",
        "G[+56]K[+42]GGK[+56]GLGK[+56]GGAK[+56]R": "P62805[aK@6](HIST1H4A)",
        "G[+56]K[+56]GGK[+56]GLGK[+42]GGAK[+56]R": "P62805[aK@13](HIST1H4A)",
        "G[+56]K[+56]GGK[+56]GLGK[+56]GGAK[+42]R": "P62805[aK@17](HIST1H4A)",
        "G[+56]K[+56]GGK[+42]GLGK[+42]GGAK[+56]R": "P62805[aK@9][aK@13](HIST1H4A)",
        "G[+56]K[+42]GGK[+42]GLGK[+56]GGAK[+56]R": "P62805[aK@6][aK@9](HIST1H4A)",
        "G[+56]K[+42]GGK[+56]GLGK[+56]GGAK[+42]R": "P62805[aK@6][aK@17](HIST1H4A)",
        "G[+56]K[+56]GGK[+56]GLGK[+42]GGAK[+42]R": "P62805[aK@13][aK@17](HIST1H4A)",
        "G[+56]K[+56]GGK[+42]GLGK[+42]GGAK[+42]R": "P62805[aK@9][aK@13][aK@17](HIST1H4A)",
        "G[+56]K[+42]GGK[+42]GLGK[+42]GGAK[+56]R": "P62805[aK@6][aK@9][aK@13](HIST1H4A)",
        "G[+56]K[+42]GGK[+42]GLGK[+56]GGAK[+42]R": "P62805[aK@6][aK@9][aK@17](HIST1H4A)",
        "G[+56]K[+42]GGK[+42]GLGK[+42]GGAK[+42]R": "P62805[aK@6][aK@9][aK@13][aK@17](HIST1H4A)",
        "G[+56]K[+56]GGK[+56]GLGK[+56]GGAK[+70]R": "P62805[meK@17](HIST1H4A)",
        "K[+112.1]VLR": "P62805[me0K@21](HIST1H4A)",
        "K[+126.1]VLR": "P62805[meK@21](HIST1H4A)",
        "K[+84.1]VLR": "P62805[me2K@21](HIST1H4A)",
        "K[+98.1]VLR": "P62805[me3K@21](HIST1H4A)",
        "D[+56]AVTYTEHAK[+56]R": "P62805-NORM(HIST1H4A)",
        "Y[+56]QK[+28]STELLIR": "P68431[me2K@57](HIST1H3A)"
    }
    meta_data = pd.read_csv('static/data/P100-all-plates-metaData.csv', sep=',', header=None)
    #print meta_data
    for i in range(0, 96):
        # print str(meta_data[6][29 + i])
        # print petide2ptmAndGene[str(meta_data[6][29 + i])]
        # print str(meta_data[11][29 + i])
        # print str(meta_data[2][29 + i])
        row_peptide.append('Peptide: ' + str(meta_data[6][29 + i]))
        shorthand = peptide2ptmAndGene[str(meta_data[6][29 + i])]
        row_shorthand.append('PTM Proteins: ' + shorthand)
        row_protein.append('Protein: ' + str(meta_data[11][29 + i]))
        row_gene.append('Gene: ' + str(meta_data[2][29 + i]))
        print i

    col_cell_line = []
    col_pert_dose = []
    col_pert_time = []
    col_pert_name = []
    col_pert_signature = []

    col_pert_type = []
    first_loop = True




    # for i in range(0, 96):
    #     row_peptide.append('Peptide: ' + str(P100_data[6][29 + i]))
    #     shorthand = petide2ptmAndGene[str(P100_data[6][29 + i])]
    #     row_shorthand.append('PTM Proteins: ' + shorthand)
    #     row_protein.append('Protein: ' + str(P100_data[11][29 + i]))
    #     row_gene.append('Gene: ' + str(P100_data[2][29 + i]))

    # open('static/data/P100_processed_perturb_for_clustergramm_slicing.json') as f:

    with open('static/data/P100_processed_perturb_for_clustergramm_slicing.json') as f:
        data = json.load(f)
        cp_iter = 1
        for var in compounds:
            cp = str(var).upper()
            print cp
            if cp in data:
                for iter in data[cp]:

                    print cp_iter
                    col_cell_line.append('Cell Line: ' + iter["Cell_line"])
                    col_pert_dose.append('Dose: ' + iter["Dose"])
                    col_pert_time.append('Time: ' + iter["Time"])
                    col_pert_name.append('Perturbations: ' + iter["Perturbations"])
                    col_pert_signature.append('Signature ID: ' + iter["Sig_id"])
                    col_pert_type.append('Type: ' + iter["Type"])

                    x1 = np.array(iter["data"])

                    # print "end of loop1"
                    x2 = x1[:, np.newaxis]
                    # print "end of loop2"
                    my_data_final = np.asarray(x2)
                    # print "end of loop3"
                    if first_loop:
                        my_data_con = my_data_final
                        first_loop = False

                    else:
                        my_data_con = np.hstack((my_data_con, my_data_final))


        # for i in range(0, 96):
        #     row_peptide.append('Peptide: ' + str(P100_data[6][29 + i]))
        #     shorthand = petide2ptmAndGene[str(P100_data[6][29 + i])]
        #     row_shorthand.append('PTM Proteins: ' + shorthand)
        #     row_protein.append('Protein: ' + str(P100_data[11][29 + i]))
        #     row_gene.append('Gene: ' + str(P100_data[2][29 + i]))

        arrays_rows2 = [np.array(row_shorthand), np.array(row_peptide), np.array(row_protein), np.array(row_gene)]
        arrays_columns2 = [np.array(col_pert_name), np.array(col_pert_signature), np.array(col_pert_time),
                           np.array(col_pert_dose), np.array(col_pert_type), np.array(col_cell_line)]
        # print arrays_columns2
        # arrays_rows2 = np.load("static/data/P100_rowMetaDat")

        tuples_rows2 = list(zip(*arrays_rows2))
        # print tuples_rows2
        tuples_columns2 = list(zip(*arrays_columns2))
        # print tuples_columns2
        rows_labels2 = pd.MultiIndex.from_tuples(tuples_rows2)
        columns_labels2 = pd.MultiIndex.from_tuples(tuples_columns2)

        # my_data_con = my_data_con.fillna(0)
        print "after na to zero"
        a2 = np.nan_to_num(my_data_con)
        # a2 = my_data_con
        aa = a2.astype(np.float)

        print "a2.shape"
        print a2.shape
        print "columns_labels2.shape"
        print columns_labels2.shape
        print "rows_labels2.shape"
        print rows_labels2.shape
        # {'nop': row1, 'o0p': row2, 'zaz': row3, 'zax': row4, 'oof': row5, 'oye': row6}
        # df = pd.DataFrame(a, index=rows_labels, columns=columns_labels)
        df2 = pd.DataFrame(aa, index=rows_labels2, columns=columns_labels2)



            #print results
    print "finished"
    P100_all_json = loadFile.make_json_from_txt(df2)
    return P100_all_json













@app.route("/api/clust/")
def make_clustergram():
    # header1 = pd.MultiIndex.from_product([['a', 'b', 'c'], ['S1', 'S2']], names=['Col', 'S'])
    # header2 = pd.MultiIndex.from_product([['a', 'b', 'c'], ['T1', 'T2']], names=['Row', 'T'])
    # df = pd.DataFrame(np.random.randint(low=1, high=10, size=(6,6)), columns = header1, index = header2)
    # df.rename_axis("Column", axis="columns")
    # df.rename_axis("Row", axis="index")

    # arrays = [['bar', 'bar', 'baz', 'baz', 'foo', 'foo', 'qux', 'qux'],
    #           ['one', 'two', 'one', 'two', 'one', 'two', 'one', 'two'],
    #           ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']]
    #columns = pd.MultiIndex.from_tuples(tuples1, names=['first', 'second'], tuples2, names=['first', 'second'])
    # s = pd.Series(np.random.randn(8), index=index)
    #df = pd.DataFrame(np.random.randint(low=1, high=10, size=(8, 8)), columns=columns, index=index)

    arrays_rows = [np.array(['Row Label: bar1', 'Row Label: arb1', 'Row Label: baz1', 'Row Label: caz1', 'Row Label: footttttt', 'Row Label: fos1']),  # Row headers
              np.array(['one', 'two', 'one', 'two', 'one', 'two', 'one', 'two']),  # Category 1
              np.array(['Letters: a', 'Letters: b', 'Letters: c', 'Letters: d', 'Letters: e', 'Letters: f', 'Letters: g', 'Letters: h']),  # Category 2
              np.array(['Numbers: 1', 'Numbers: 2', 'Numbers: 3', 'Numbers: 1', 'Numbers: 2', 'Numbers: 3', 'Numbers: 1', 'Numbers: 2'])]  # Category 3
    arrays_columns = [np.array(['Column Label: nop', 'Column Label: o0p', 'Column Label: zaz', 'Column Label: zax', 'Column Label: oof', 'Column Label: oye']),  # Column headers
              np.array(['Symbol: @', 'Symbol: !', 'Symbol: #', 'Symbol: @', 'Symbol: !', 'Symbol: #', 'Symbol: @', 'Symbol: !']),  # Category 1
              np.array(['three', 'four', 'three', 'four', 'three', 'four', 'three', 'four'])]  # Category 2
    # s = pd.Series(np.random.randn(8), index=arrays)
    # df = pd.DataFrame(np.random.randn(3, 8), index=['A', 'B', 'C'], columns=index)
    tuples_rows = list(zip(*arrays_rows))
    tuples_columns = list(zip(*arrays_columns))
    rows_labels = pd.MultiIndex.from_tuples(tuples_rows)
    columns_labels = pd.MultiIndex.from_tuples(tuples_columns)

    a = np.matrix('.1 .2 .3 .4 .5 .6;'
                  ' .2 .3 .4 .5 .6 .7;'
                  ' .3 .4 .5 .6 .7 .8;'
                  ' -.1, -.2 -.3 -.4, -.5, -.6;'
                  ' -.2, -.3 -.4 -.5, -.6, -.7;'
                  ' -.3, -.4 -.5 -.6, -.7, -.8')

    # {'nop': row1, 'o0p': row2, 'zaz': row3, 'zax': row4, 'oof': row5, 'oye': row6}
    df = pd.DataFrame(a, index=rows_labels, columns=columns_labels)
    # df.rename_axis("Row", axis="index")

    # df = pd.DataFrame(np.random.randn(6, 6), index=rows_labels[:6], columns=columns_labels[:6])
    #with open("main/Service/rc_two_cats.txt") as myfile:
        # data = myfile.read()
    #    data = pd.read_csv(myfile)

    # data = json.loads(response.text)
    # print(data)
    # string = df.astype(str)
    # return string
    # return loadFile.make_json_from_txt(data)
    example_json = loadFile.make_json_from_txt(df)
    json.dump(example_json, open("static/data/example_clustergram.json", 'w'))
    return example_json


@app.route("/version")
def version():
    return make_response(open('static/versions-mapping.json').read())

#@app.route("/hello")
#def hello():
#    return "Hello World!"

@app.route("/api/submitGeneToIlincs/<gene_name>")# example: bcl2a1
def submitGeneToIlincsMethod(gene_name):
    data = ilincsSearch.search_gene(gene_name)
    #data = json.load(urllib2.urlopen(link))
    #return str(data)
    sigid = []

    for item in data:
        sigid.append(str(item['signatureid']))


    # printing all signatures associated with gene
    print sigid

    # getting concordant signatures for each id
    for id in sigid:
        con = []
        dis = []
        #link = "http://www.ilincs.org/api/SignatureMeta/findConcordantSignatures?sigID=%22" + str(id) + "%22&lib=%22LIB_5%22"
        #data = json.load(urllib2.urlopen(link))
        data = ilincsSearch.search_sigid(id)
        for perturbagen in data:
            pid = str(perturbagen['perturbagenID'])
            sim = str(perturbagen['similarity'])

            # concordant
            if sim > 0.321:
                con.append(pid)
            # discordant
            elif sim < -0.235:
                dis.append(pid)
    return str(con)

        #print "\nsignature ---- > ", id, "concordant ----->", con, "discordant ------>", dis




@app.route('/favicon.ico')
def favicon():
    return send_from_directory(os.path.join(app.root_path, 'static'),
                               'img/favicon.ico')


# if __name__ == '__main__':
#     app.run(debug=True)





if __name__ == "__main__":
    app.run()
