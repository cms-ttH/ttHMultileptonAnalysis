#!/bin/bash

echo ""
echo "-----------------"
echo "Generating limits"
echo "-----------------"
echo ""


## Blind, real data
makeCardParams='--jet_tag --label'
#makeCardParams='--jet_tag --label'

#combinePrefix='-M HybridNew --testStat=PL --singlePoint 0 --onlyTestStat'
## Bkg+SM sig, use data: --toysFreq
## Verbose (pulls): -v 2
#combinePrefix='-M HybridNew --testStat=PL --singlePoint 0 --onlyTestStat --toysFreq -v 2'
#combinePrefix='-M MaxLikelihoodFit --minos all --toysFreq'
combinePrefix='-M MaxLikelihoodFit --minos all'

combineSuffix='-t -1 --expectSignal=1'
#combineSuffix='-t -1 --expectSignal=1'

# ## Unblind, fake data
# makeCardParams='--data_sig_bkg --jet_tag --label'
# combinePrefix='-M HybridNew --testStat=PL --singlePoint 0 --onlyTestStat --toysFreq'
# combineSuffix='-m 125'

# ## October 31 parameters
# systFile=systsEnv_Sep10.txt
# label_=Sep10_
# systFile=systsEnv_Oct31.txt
# label_=Oct31_no_bin_stat_

## November 14 parameters
systFile=systsEnv_Nov14.txt
label_=Nov14_some_bin_stat_

SSttWvar=FinalBDT
threeLttWvar=FinalBDT
ttWvar=FinalBDT
threeLttZvar=FinalBDT
fourLttZvar=numMediumBJets
OSttZvar=BDT_zjets_ttbar

# SSttWvar=FinalBDT_Oct31
# threeLttWvar=FinalBDT_Oct31
# ttWvar=FinalBDT_Oct31
# threeLttZvar=FinalBDT_Oct31
# fourLttZvar=numMediumBJets
# OSttZvar=BDT_zjets_ttbar_Oct31


echo ""
echo "-----------------"
echo "Making shape cards"
echo "-----------------"
echo ""

python makeShapeCards.py yaml/limit_configuration_SS_ttW_2t_lepCut.yaml $SSttWvar $systFile $makeCardParams 2t_$label_$SSttWvar
python makeShapeCards.py yaml/limit_configuration_3l_ttW_2t_lepCut.yaml $threeLttWvar $systFile $makeCardParams 2t_$label_$threeLttWvar
python makeShapeCards.py yaml/limit_configuration_3l_ttZ_2t_lepCut.yaml $threeLttZvar $systFile $makeCardParams 2t_$label_$threeLttZvar
#python makeShapeCards.py yaml/limit_configuration_4l_ttZ_2t_lepCut.yaml $fourLttZvar $systFile $makeCardParams 2t_$label_$fourLttZvar
python makeShapeCards.py yaml/limit_configuration_4l_ttZ_4l_lepCut.yaml $fourLttZvar $systFile $makeCardParams 4l_$label_$fourLttZvar
python makeShapeCards.py yaml/limit_configuration_OS_ttZ_2l_lepCut.yaml $OSttZvar $systFile $makeCardParams 2l_$label_$OSttZvar

echo ""
echo "-----------------"
echo "Combining cards"
echo "-----------------"
echo ""

combineCards.py ND_cards/ttW_2lss_*_*j_bloose_2t_$label_$SSttWvar.card.txt > ND_cards/cardComb/ttW_2lss_bloose_2t_$label_$SSttWvar.card.txt
combineCards.py ND_cards/ttW_3l_*j_bloose_2t_$label_$threeLttWvar.card.txt > ND_cards/cardComb/ttW_3l_bloose_2t_$label_$threeLttWvar.card.txt
combineCards.py ND_cards/ttZ_3l_*j_bloose_2t_$label_$threeLttZvar.card.txt > ND_cards/cardComb/ttZ_3l_bloose_2t_$label_$threeLttZvar.card.txt
#combineCards.py ND_cards/ttZ_4l_ge1j_Z*_mht30_1bloose_2t_$label_$fourLttZvar.card.txt > ND_cards/cardComb/ttZ_4l_ge1j_mht30_1bloose_2t_$label_$fourLttZvar.card.txt
combineCards.py ND_cards/ttZ_4l_ge1j_Z*_mht30_1bloose_4l_$label_$fourLttZvar.card.txt > ND_cards/cardComb/ttZ_4l_ge1j_mht30_1bloose_4l_$label_$fourLttZvar.card.txt
combineCards.py ND_cards/ttZ_2los_*j_*t_2l_$label_$OSttZvar.card.txt > ND_cards/cardComb/ttZ_2los_2l_$label_$OSttZvar.card.txt

combineCards.py ND_cards/ttW_2lss_*_*j_bloose_2t_$label_$SSttWvar.card.txt ND_cards/ttW_3l_*j_bloose_2t_$label_$threeLttWvar.card.txt > ND_cards/cardComb/ttW_23l_2t_$label_$ttWvar.card.txt
#combineCards.py ND_cards/ttZ_3l_*j_bloose_2t_$label_$threeLttZvar.card.txt ND_cards/ttZ_4l_ge1j_Z*_mht30_1bloose_2t_$label_$fourLttZvar.card.txt > ND_cards/cardComb/ttZ_34l_2t_$label_$threeLttZvar.card.txt
#combineCards.py ND_cards/ttZ_3l_*j_bloose_2t_$label_$threeLttZvar.card.txt ND_cards/ttZ_4l_ge1j_Z*_mht30_1bloose_2t_$label_$fourLttZvar.card.txt ND_cards/ttZ_2los_*j_*t_2l_$label_$OSttZvar.card.txt > ND_cards/cardComb/ttZ_234l_$label_$threeLttZvar.card.txt
combineCards.py ND_cards/ttZ_3l_*j_bloose_2t_$label_$threeLttZvar.card.txt ND_cards/ttZ_4l_ge1j_Z*_mht30_1bloose_4l_$label_$fourLttZvar.card.txt > ND_cards/cardComb/ttZ_34l_$label_$threeLttZvar.card.txt
combineCards.py ND_cards/ttZ_3l_*j_bloose_2t_$label_$threeLttZvar.card.txt ND_cards/ttZ_4l_ge1j_Z*_mht30_1bloose_4l_$label_$fourLttZvar.card.txt ND_cards/ttZ_2los_*j_*t_2l_$label_$OSttZvar.card.txt > ND_cards/cardComb/ttZ_234l_$label_$threeLttZvar.card.txt

# ## Apparently this isn't how you combine two different signals
# ##combineCards.py ND_cards/ttW_2lss_*_*j_bloose_2t_$label_$SSttWvar.card.txt ND_cards/ttW_3l_*j_bloose_2t_$label_$threeLttWvar.card.txt  ND_cards/ttZ_3l_*j_bloose_2t_$label_$threeLttZvar.card.txt ND_cards/ttZ_4l_ge1j_Z*_mht30_1bloose_2t_$label_$fourLttZvar.card.txt > ND_cards/cardComb/ttW_23l_ttZ_34l_2t_$label_.card.txt

echo ""
echo "-----------------"
echo "Computing limits"
echo "-----------------"
echo ""

## get proper version of combine
cd /afs/crc.nd.edu/user/a/abrinke1/CMSSW_6_1_1/src/HiggsAnalysis/CombinedLimit
eval `scramv1 runtime -sh`
cd -

echo "-----------------"
echo "SS ttW limits - "$SSttWvar
echo "-----------------"
combine $combinePrefix ND_cards/cardComb/ttW_2lss_bloose_2t_$label_$SSttWvar.card.txt $combineSuffix
echo "-----------------"
echo "3l ttW limits - "$threeLttWvar
echo "-----------------"
combine $combinePrefix ND_cards/cardComb/ttW_3l_bloose_2t_$label_$threeLttWvar.card.txt $combineSuffix
echo "-----------------"
echo "3l ttZ limits - "$threeLttZvar
echo "-----------------"
combine $combinePrefix ND_cards/cardComb/ttZ_3l_bloose_2t_$label_$threeLttZvar.card.txt $combineSuffix
echo "-----------------"
echo "4l ttZ limits - "$fourLttZvar
echo "-----------------"
#combine $combinePrefix ND_cards/cardComb/ttZ_4l_ge1j_mht30_1bloose_2t_$label_$fourLttZvar.card.txt $combineSuffix
combine $combinePrefix ND_cards/cardComb/ttZ_4l_ge1j_mht30_1bloose_4l_$label_$fourLttZvar.card.txt $combineSuffix
echo "-----------------"
echo "OS ttZ limits - "$OSttZvar
echo "-----------------"
combine $combinePrefix ND_cards/cardComb/ttZ_2los_2l_$label_$OSttZvar.card.txt $combineSuffix
echo "-----------------"
echo "SS and 3l ttW limits - "$ttWvar
echo "-----------------"
combine $combinePrefix ND_cards/cardComb/ttW_23l_2t_$label_$ttWvar.card.txt $combineSuffix
echo "-----------------"
echo "3l and 4l ttZ limits - "$threeLttZvar
echo "-----------------"
#combine $combinePrefix ND_cards/cardComb/ttZ_34l_2t_$label_$threeLttZvar.card.txt $combineSuffix
combine $combinePrefix ND_cards/cardComb/ttZ_34l_$label_$threeLttZvar.card.txt $combineSuffix
echo "-----------------"
echo "3l, 4l, and OS ttZ limits - "$threeLttZvar
echo "-----------------"
combine $combinePrefix ND_cards/cardComb/ttZ_234l_$label_$threeLttZvar.card.txt $combineSuffix

eval `scramv1 runtime -sh`

# grep -H -r "ln Q_{Profile}" log.txt | cut -f2- -d'='
# echo 'scale=30;sqrt(3)' | bc
