#!/usr/bin/env python
import os, sys
from argparse import ArgumentParser
import copy
import math
import yaml
import ttHMultileptonAnalysis.DrawPlots.utilities.plot_helper as plot_helper
from ttHMultileptonAnalysis.DrawPlots.utilities.yamlhelpers import ttHMultileptonYAMLLoader
from ttHMultileptonAnalysis.DrawPlots.utilities.prettytable import PrettyTable
from ttHMultileptonAnalysis.DrawPlots.utilities.ordereddict import DefaultOrderedDict
from ttHMultileptonAnalysis.DrawPlots.utilities.pie import make_yield_pie_charts

parser = ArgumentParser(description='Make stack plots from histogram files.')
parser.add_argument('config_file_name', nargs='?', default='stack_plot_configuration.yaml', help='Configuration file to process.')
parser.add_argument('cosmetics_file_name', nargs='?', default='stack_plot_cosmetics.yaml', help='Cosmetics configuration file to process.')
parser.add_argument('-w', '--web', action='store_true', help='post each plot to the user\'s AFS space')
parser.add_argument('--label', help='Override the label designated in the configuration file with LABEL')
parser.add_argument('-p', '--pie', action='store_true', help='Post pie plots of background composition and print signal/bg table')
parser.add_argument('--blind', action='store_true', help='Make data-blind plots')
parser.add_argument('--postFit', action='store_true', help='Use post-fit histograms')
parser.add_argument('--postFitNorms', action='store_true', help='Use post-fit normalizations and uncertainties, but not shapes')
args = parser.parse_args()

## config is mostly options that change the substance of what is displayed
## config_cosmetics is mostly options that change the style of the stack plot
with open(args.config_file_name) as config_file:
    config = yaml.load(config_file, Loader=ttHMultileptonYAMLLoader)

with open(args.cosmetics_file_name) as cosmetics_file:
    cosmetics = yaml.load(cosmetics_file, Loader=ttHMultileptonYAMLLoader)

lepton_categories = config['lepton categories']
jet_tag_categories = config['jet tag categories']
distributions = config['distributions']
if args.postFit:
    distributions = config['postFit distributions']
if args.postFitNorms:
    for distribution in config['postFit distributions']:
        distribution_norm =  distribution
#distributions['integral_histo'] = ['isCleanEvent', False, False]
signal_samples = config.get('signal samples', {})
background_samples = config.get('background samples', {})
draw_names = lambda d: [d[x]['draw name'] for x in d.keys()]
yields = plot_helper.Yields(jet_tag_categories, draw_names(signal_samples)+draw_names(background_samples))
raw_yields = plot_helper.Yields(jet_tag_categories, draw_names(signal_samples)+draw_names(background_samples))

def main():
    ## import the root libraries; with this import you won't need to prefix your root objects with ROOT
    from ROOT import *
    gROOT.SetBatch()

    if args.label:
        config['input file label'] = args.label

    out_location = config['output file location']
    if args.blind:
        out_location += '_blind'
    if args.postFit:
        out_location += '_postFitS'
    if args.postFitNorms:
        out_location += '_postFitSNorms'
    if args.cosmetics_file_name.startswith('cosmetics_'):
        out_location += (args.cosmetics_file_name.replace('cosmetics','')).replace('.yaml','')

    if args.web:
        www_plot_directories = []
        for lepton_category in lepton_categories:
            for jet_tag_category in jet_tag_categories:
                www_plot_directories.append(os.path.join(out_location, lepton_category, jet_tag_category))

        plot_helper.setup_web_posting(www_plot_directories, 4, args.config_file_name, args.cosmetics_file_name)

    for lepton_category in lepton_categories:
        print '\n\nStarting lepton category %s...\n' % lepton_category

        for jet_tag_category in jet_tag_categories:
            print 'Starting jet tag category %s...' % jet_tag_category

            for distribution in distributions:
                print 'Drawing distribution: %s with jet selection printing name: %s' %  (distribution, jet_tag_categories[jet_tag_category])
                draw_stack_plot(lepton_category, jet_tag_category, distribution)

            #draw_stack_plot(lepton_category, jet_tag_category, 'integral_histo')

    print_yield_table(yields, 'final yields', '.2f')
    print_yield_table(raw_yields, 'raw yields', '.0f')

    if args.pie:
        make_yield_pie_charts(yields, out_location, draw_names(signal_samples), jet_tag_categories)

    if args.web:
        plot_helper.update_indexes(os.path.join(out_location, '..'))
        print '\nFinished processing.  Plots will be posted to: http://www.crc.nd.edu/~%s/%s' % (os.environ['USER'], out_location)

def print_yield_table(yields, title, precision):
    columns = ['sample'] + lepton_categories.keys()
    for jet_tag_category in jet_tag_categories:
        print "\n\n%s: (jet tag category %s)" % (title, jet_tag_categories[jet_tag_category])
        yield_table = PrettyTable(columns)
        for sample in yields[jet_tag_category]:
            entry = [sample]
            for lepton_category in lepton_categories:
                if yields[jet_tag_category][sample].has_key(lepton_category):
                    entry.append(format(yields[jet_tag_category][sample][lepton_category], precision))
                else:
                    entry.append('-')
            yield_table.add_row(entry)
        print yield_table

def draw_stack_plot(lepton_category, jet_tag_category, distribution):
    stack_plot = THStack("theStack", "")
    stack_plot_legend = make_legend()

    (luminosity_info_tex, selection_info_tex, CMS_info_tex, fit_status_tex, SF_info_tex) = make_info_tex_objects(lepton_category, jet_tag_category)

    out_location = config['output file location']
    if args.blind:
        out_location = config['output file location']+'_blind'
    if args.postFit:
        out_location = config['output file location']+'_postFitS'
    if args.postFitNorms:
        out_location = config['output file location']+'_postFitSNorms'
    if args.cosmetics_file_name.startswith('cosmetics_'):
        out_location += (args.cosmetics_file_name.replace('cosmetics','')).replace('.yaml','')
        
    histogram_dictionary = {}
    mc_sum = 0.0
    signal_sum = 0.0
    data_sum = 0.0
    all_signals = 0.0

    systematics_by_sample = {}
    # Sum histograms in each sample group, then add summed histos to stack plot, legend, and histogram_dictionary
    for sample_group in background_samples: # Add samples in each sample group together
        if sample_group in lepton_categories[lepton_category].get('excluded samples', []):
            continue
        if 'sideband' in sample_group and not plot_helper.is_matching_data_sample(lepton_categories[lepton_category].get('data samples', []), sample_group):
            continue
        systematics = plot_helper.customize_systematics(config['systematics'], background_samples[sample_group].get('systematics', 'common'))
        if config['skip systematics']:
            systematics = ['nominal']
        systematics_by_sample[sample_group] = copy.deepcopy(systematics)
        for systematic in systematics:
            samples_in_group = background_samples[sample_group]['samples']
            if args.postFit:
               samples_in_group = background_samples[sample_group]['postFit samples'] 
            group_histogram = get_group_histogram(distribution, systematic, samples_in_group, lepton_category, jet_tag_category, sample_group)
            # Add one entry (the sum of all samples in the group) for each sample group to the
            # dictionary; the group will be treated as a single sample for subsequent error calculations
            if args.postFitNorms and (group_histogram or 'nominal' in systematic):
               samples_in_group_norms = background_samples[sample_group]['postFit samples']
               args.postFit = true
               group_histogram_norm = get_group_histogram(distribution_norm, 'nominal', samples_in_group_norms, lepton_category, jet_tag_category, sample_group)                                              

               args.postFit = false
               if group_histogram_norm.Integral() > 0 and group_histogram.Integral() > 0:
                   group_histogram.Scale(group_histogram_norm.Integral()/group_histogram.Integral())

                   if 'nominal' in systematic:
                       relative_error = group_histogram_norm.GetBinError(1) / group_histogram_norm.GetBinContent(1)
                       for bin in range(1, group_histogram.GetNbinsX()+1):
                           if group_histogram.GetBinContent(bin) > 0:
                               old_bin_error = group_histogram.GetBinError(bin)
                               group_histogram.SetBinError(bin, math.sqrt( pow(old_bin_error, 2) + pow(relative_error*group_histogram.GetBinContent(bin), 2) ))                   
               else:
                   if (group_histogram_norm.Integral() > 0.01) + (group_histogram.Integral() > 0.01) != 0:
                       print 'ERROR! %s has pre-fit value of %.2f, post-fit of %.2f' % (sample_group, group_histogram.Integral(), group_histogram_norm.Integral())
               
            try:
                histogram_dictionary[sample_group+'_'+systematic] = group_histogram.Clone()
            except AttributeError:
                print 'Problems finding any input files for sample group %s, systematic %s.  Skipping sample group / systematic ...' % (sample_group, systematic)
                if systematic in systematics_by_sample[sample_group]:
                    systematics_by_sample[sample_group].remove(systematic)
                continue
            if systematic == 'nominal':
                histogram = group_histogram.Clone('stack')
                mc_sum += histogram.Integral()
                stack_plot.Add(histogram, 'hist')
                if cosmetics['legend precision'] == -1:
                    stack_plot_legend.AddEntry(histogram, '%s' % (background_samples[sample_group]['draw name']), 'f' )
                elif cosmetics['legend precision'] == 0 or histogram.Integral() >= 19.95:
                    stack_plot_legend.AddEntry(histogram, '%s (%.0f)' % (background_samples[sample_group]['draw name'], histogram.Integral()), 'f')
                else:
                    stack_plot_legend.AddEntry(histogram, '%s (%0.1f)' % (background_samples[sample_group]['draw name'], histogram.Integral()), 'f')

                yields[jet_tag_category][background_samples[sample_group]['draw name']][lepton_category] = histogram.Integral()
                raw_yields[jet_tag_category][background_samples[sample_group]['draw name']][lepton_category] = histogram.GetEntries()

    ## Draw the signal sample histogram(s), put in legend
    signal_histograms = {}
    for sample_group in signal_samples:
        systematics = plot_helper.customize_systematics(config['systematics'], signal_samples[sample_group].get('systematics', 'common'))
        if config['skip systematics']:
            systematics = ['nominal']
        systematics_by_sample[sample_group] = copy.deepcopy(systematics)
        for systematic in systematics:
            samples_in_group = signal_samples[sample_group]['samples']
            if args.postFit:
                samples_in_group = signal_samples[sample_group]['postFit samples']
            group_histogram = get_group_histogram(distribution, systematic, samples_in_group, lepton_category, jet_tag_category, sample_group)
            if args.postFitNorms and (group_histogram or 'nominal' in systematic):
                samples_in_group_norms = signal_samples[sample_group]['postFit samples']
                args.postFit = true
                group_histogram_norm = get_group_histogram(distribution_norm, 'nominal', samples_in_group_norms, lepton_category, jet_tag_category, sample_group)
 
                args.postFit = false
                if group_histogram_norm.Integral() > 0 and group_histogram.Integral() > 0:
                    group_histogram.Scale(group_histogram_norm.Integral()/group_histogram.Integral())

                    if 'nominal' in systematic:
                        relative_error = group_histogram_norm.GetBinError(1) / group_histogram_norm.GetBinContent(1)
                        for bin in range(1, group_histogram.GetNbinsX()+1):
                            if group_histogram.GetBinContent(bin) > 0:
                                old_bin_error = group_histogram.GetBinError(bin)
                                group_histogram.SetBinError(bin, math.sqrt( pow(old_bin_error, 2) + pow(relative_error*group_histogram.GetBinContent(bin), 2) ))                   
                    
                else:
                    if (group_histogram_norm.Integral() > 0.01) + (group_histogram.Integral() > 0.01) != 0:
                        print 'ERROR! %s has pre-fit value of %.2f, post-fit of %.2f' % (sample_group, group_histogram.Integral(), group_histogram_norm.Integral())

            try:
                histogram_dictionary[sample_group+'_'+systematic] = group_histogram.Clone()
            except AttributeError:
                print 'Problems finding any input files for sample group %s, systematic %s.  Skipping sample group / systematic ...' % (sample_group, systematic)
                if systematic in systematics_by_sample[sample_group]:
                    systematics_by_sample[sample_group].remove(systematic)
                continue
            if systematic == 'nominal':
                signal_histograms[sample_group] = group_histogram.Clone(sample_group)
                legend_option = 'l'
                if signal_samples[sample_group]['stack or line'] == 'stack':
                    stack_plot.Add(signal_histograms[sample_group], 'hist')
                    legend_option = 'f'
                elif signal_samples[sample_group]['stack or line'] != 'line':
                    sys.exit('For signal, must be line or stack')
#                 elif signal_samples[sample_group]['stack or line'] == 'line':
#                     stack_plot.Add(signal_histogram, 'hist')

                signal_sum = signal_histograms[sample_group].Integral()
                signal_scaling = 1.0
                if signal_sum > 0:
                    if signal_samples[sample_group]['scale'] == 'norm':
                        signal_scaling = mc_sum / signal_sum
                    else:
                        signal_scaling = signal_samples[sample_group]['scale']
                        
                    signal_histograms[sample_group].Scale(signal_scaling)
                    if signal_scaling == 1.0:
                        if cosmetics['legend precision'] == -1:
                            stack_plot_legend.AddEntry(signal_histograms[sample_group], '%s' % (signal_samples[sample_group]['draw name']), legend_option)
                        elif cosmetics['legend precision'] == 0 or signal_sum >= 19.95:
                            stack_plot_legend.AddEntry(signal_histograms[sample_group], '%s (%.0f)' % (signal_samples[sample_group]['draw name'], signal_sum), legend_option)
                        else:
                            stack_plot_legend.AddEntry(signal_histograms[sample_group], '%s (%0.1f)' % (signal_samples[sample_group]['draw name'], signal_sum), legend_option)
                    elif cosmetics['legend precision'] == -1:
                        if signal_scaling < 19.95:
                            stack_plot_legend.AddEntry(signal_histograms[sample_group], '%s x %0.1f' % (signal_samples[sample_group]['draw name'], signal_scaling), legend_option)
                        else:
                            stack_plot_legend.AddEntry(signal_histograms[sample_group], '%s x %.0f' % (signal_samples[sample_group]['draw name'], signal_scaling), legend_option)
                    elif cosmetics['legend precision'] == 0 or (signal_sum >= 19.95 and signal_scaling >= 19.95):
                        stack_plot_legend.AddEntry(signal_histograms[sample_group], '%s (%.0f x %.0f)' % (signal_samples[sample_group]['draw name'], signal_sum, signal_scaling), legend_option)
                    elif signal_sum >= 19.95:
                        stack_plot_legend.AddEntry(signal_histograms[sample_group], '%s (%.0f x %0.1f)' % (signal_samples[sample_group]['draw name'], signal_sum, signal_scaling), legend_option)
                    elif signal_scaling >= 19.95:
                        stack_plot_legend.AddEntry(signal_histograms[sample_group], '%s (%0.1f x %.0f)' % (signal_samples[sample_group]['draw name'], signal_sum, signal_scaling), legend_option)
                    else:
                        stack_plot_legend.AddEntry(signal_histograms[sample_group], '%s (%0.1f x %0.1f)' % (signal_samples[sample_group]['draw name'], signal_sum, signal_scaling), legend_option)

                else:
                    stack_plot_legend.AddEntry(signal_histograms[sample_group], '%s (0.0)' % signal_samples[sample_group]['draw name'], legend_option)

                yields[jet_tag_category][signal_samples[sample_group]['draw name']][lepton_category] = histogram_dictionary[sample_group+'_nominal'].Integral()
                raw_yields[jet_tag_category][signal_samples[sample_group]['draw name']][lepton_category] = histogram_dictionary[sample_group+'_nominal'].GetEntries()
                all_signals += signal_sum

    ## Draw the data histogram, put in legend
    skip_data = args.blind
    if not any([lepton_categories[lep_cat].get('data samples') for lep_cat in lepton_categories]):
        skip_data = True
    if not skip_data:
        group_histogram = None
        samples_in_group = lepton_categories[lepton_category].get('data samples', [])
        if args.postFit:
            samples_in_group = lepton_categories[lepton_category].get('postFit data samples', [])
        if args.postFitNorms:
            samples_in_group_norms = lepton_categories[lepton_category].get('postFit data samples', [])
        
        group_histogram = get_group_histogram(distribution, 'nominal', samples_in_group, lepton_category, jet_tag_category, 'data')
        if not group_histogram:
            print 'Problem finding the requested data histograms.  Skipping them and continuing on...'
            skip_data = True
        else:
            histogram_dictionary['data'] = group_histogram.Clone()
            data_histogram = group_histogram.Clone('data')
            data_sum = data_histogram.Integral()
            if cosmetics['legend precision'] == -1:
                stack_plot_legend.AddEntry(data_histogram, 'Data', 'lp')
            else:
                stack_plot_legend.AddEntry(data_histogram, 'Data (%.0f)' % data_sum, 'lp')
            yields[jet_tag_category]['data'][lepton_category] = data_histogram.Integral()
            raw_yields[jet_tag_category]['data'][lepton_category] = data_histogram.GetEntries()

    lumi_error = config['luminosity error']
    trigger_SF_error = config['trigger SF error']
    lumi_trigger_SF_error = math.sqrt(math.pow(lumi_error, 2) + math.pow(trigger_SF_error, 2))

    # Create a histogram with the error bars for the MC stack
    draw_mc_error_histo = config.get('draw mc error histo', True)
    try:
        nBins = stack_plot.GetStack().Last().GetNbinsX()
        xMin = stack_plot.GetStack().Last().GetXaxis().GetXmin()
        xMax = stack_plot.GetStack().Last().GetXaxis().GetXmax()
    except ReferenceError:
        if len(histogram_dictionary.keys()) > 0:
            draw_mc_error_histo = False
        else:
            print 'No histograms in stack for distribution %s.  Skipping it and continuing on...' % distribution
            return

    if draw_mc_error_histo:
        mc_error_histogram = stack_plot.GetStack().Last().Clone('mc_error_histogram')
        mc_error_histogram.SetTitle('')
        for i in range(1, nBins+1):
            mc_error_histogram.SetBinContent(i, stack_plot.GetStack().Last().GetBinContent(i))
            bin_error_squared = math.pow(mc_error_histogram.GetBinError(i), 2) ## Bin statistical error and cross section error
            bin_error_squared += math.pow(lumi_trigger_SF_error * stack_plot.GetStack().Last().GetBinContent(i), 2)
            ## Correlate systematics
            bin_error_syst = {}
            for sample_group, sample_systs in systematics_by_sample.items(): #systematics_by_sample is a dictionary (keys: samples and sample groups, values: systematics list)
                if sample_group in lepton_categories[lepton_category].get('excluded samples', []):
                    continue
                if 'sideband' in sample_group and not plot_helper.is_matching_data_sample(lepton_categories[lepton_category].get('data samples', []), sample_group):
                    continue
                for syst in sample_systs:
                    if not syst in bin_error_syst.keys():
                        bin_error_syst[syst] = histogram_dictionary[sample_group+'_'+syst].GetBinContent(i) - histogram_dictionary[sample_group+'_nominal'].GetBinContent(i)
                    else:
                        bin_error_syst[syst] += histogram_dictionary[sample_group+'_'+syst].GetBinContent(i) - histogram_dictionary[sample_group+'_nominal'].GetBinContent(i)

                    ## Acting as if all systematics were uncorrelated
                    #bin_error_squared += math.pow(histogram_dictionary[sample_group+'_'+systematic].GetBinContent(i) - histogram_dictionary[sample_group+'_nominal'].GetBinContent(i), 2)

            for syst in bin_error_syst.keys():
                bin_error_squared += math.pow(bin_error_syst[syst], 2)
            mc_error_histogram.SetBinError(i, math.sqrt(bin_error_squared))
            
        mc_error_histogram.SetFillStyle(cosmetics['mc error fill style'])
        mc_error_histogram.SetFillColor(cosmetics['mc error fill color'])

        if cosmetics['legend precision'] == -1:
            True
        elif cosmetics['legend precision'] == 0 or mc_sum + signal_sum >= 19.95:
            stack_plot_legend.AddEntry(mc_error_histogram, 'Sum Pred. (%.0f)' % (mc_sum + signal_sum), 'f')
        else:
            stack_plot_legend.AddEntry(mc_error_histogram, 'Sum Pred. (%0.1f)' % (mc_sum + signal_sum), 'f')

    if args.postFit:
        ## Post-fit error histogram has correct error bands but incorrect central value; need to fill legend before
        mc_error_histogram = get_histogram(distribution, 'nominal', 'err_hist_1_sigma', lepton_category, jet_tag_category)
        mc_error_histogram.SetFillStyle(cosmetics['mc error fill style'])
        mc_error_histogram.SetFillColor(cosmetics['mc error fill color'])
        mc_error_histogram_2_sigma = get_histogram(distribution, 'nominal', 'err_hist_2_sigma', lepton_category, jet_tag_category)

    yields[jet_tag_category]['all backgrounds'][lepton_category] = mc_sum
    yields[jet_tag_category]['all signals'][lepton_category] = all_signals

    plot_max = stack_plot.GetMaximum()
    if len(signal_histograms.keys()) > 0:
        plot_max = max(plot_max, max([hist.GetMaximum() for hist in signal_histograms.values()]))
    if not skip_data:
        plot_max = max( plot_max, data_histogram.GetMaximum() + math.sqrt(data_histogram.GetMaximum()) )

    stack_plot = configure_stack(stack_plot, plot_max)
    canvas = TCanvas(distribution+'Lin', distribution, cosmetics['canvas min'], cosmetics['canvas max'])

    gStyle.SetPadBorderMode(cosmetics['pad border mode'])
    gStyle.SetFrameBorderMode(cosmetics['frame border mode'])

    top_canvas, bottom_canvas = get_configured_canvases()
    top_canvas.Draw()
    bottom_canvas.Draw()

    top_canvas.cd()
    if config['log scale']:
        gPad.SetLogy()

    gPad.SetBottomMargin(cosmetics['pad bottom margin'])

    stack_plot.Draw()

    ## For some reason can't be done before Draw()
    stack_plot.GetYaxis().SetTitleSize(cosmetics['stack y axis title size'])
    stack_plot.GetYaxis().SetTitleOffset(cosmetics['stack y axis title offset'])
    stack_plot.GetYaxis().SetNdivisions(cosmetics['stack y axis Ndivisions'])
    stack_plot.GetYaxis().SetLabelSize(cosmetics['stack y axis label size'])
    stack_plot.GetYaxis().SetLabelOffset(cosmetics['stack y axis label offset'])

    stack_plot.GetXaxis().SetNdivisions(cosmetics['ratio hist x axis Ndivisions'])    
    
    if draw_mc_error_histo or args.postFit or args.postFitNorms:
        mc_error_histogram.Draw(cosmetics['mc error histogram draw style'])
        if args.postFit and abs(1 - mc_error_histogram.Integral()/stack_plot.GetStack().Last().Integral()) > 0.01:
            print 'WARNING: mc_error and stack differ by %0.2f percent overall' % (abs(1 - mc_error_histogram.Integral()/stack_plot.GetStack().Last().Integral())*100.0)
    for sample_group in signal_samples:
        if config['signal samples'][sample_group]['stack or line'] == 'line':
            signal_histograms[sample_group].Draw(cosmetics['signal histogram draw style'])

    ## asymmetrical poisson errors for data
    if not skip_data:
        ggg = get_configured_data_asymmetric_errors(data_histogram)
        ggg.Draw(cosmetics['ggg draw style'])

    ## calculate the KS test result, put it somewhere
    if config['KS test'] and not skip_data and (draw_mc_error_histo or args.postFit or args.postFitNorms):
        ks_result = data_histogram.KolmogorovTest(mc_error_histogram)
        luminosity_info_tex.SetTitle('%s (KS = %0.2f)' % (luminosity_info_tex.GetTitle(), ks_result))

    if config['draw legend']:
        stack_plot_legend.Draw()

    luminosity_info_tex.DrawLatex(cosmetics['lumi text x'], cosmetics['lumi text y'], luminosity_info_tex.GetTitle())
    if config['selection info']:
        selection_info_tex.DrawLatex(cosmetics['selection text x'], cosmetics['selection text y'], selection_info_tex.GetTitle())
        CMS_info_tex.DrawLatex(cosmetics['CMS text x'], cosmetics['CMS text y'], CMS_info_tex.GetTitle())
        fit_status_tex.DrawLatex(cosmetics['fit status text x'], cosmetics['fit status text y'], fit_status_tex.GetTitle())

    if config['SF info']:
        if (dist.find('CFMlpANN') != -1):
            SF_info_tex.DrawLatex(cosmetics['SF info min'], cosmetics['SF info max'], SF_info_tex.GetTitle())

    bottom_canvas.cd()

    if (draw_mc_error_histo or args.postFit or args.postFitNorms) and config.get('draw ratio plot', True):
        if not skip_data:
            (ratio_histogram, ratio_error_histogram) = make_ratio_histogram(nBins, xMin, xMax, mc_error_histogram, data_histogram, stack_plot)
        else:
            (ratio_histogram, ratio_error_histogram) = make_ratio_histogram(nBins, xMin, xMax, mc_error_histogram, signal_histograms.items()[0][1], stack_plot)
            ratio_histogram.Scale(1.0/signal_scaling)

        ratio_histogram = configure_ratio_histogram(ratio_histogram, distribution)
        ratio_error_histogram = configure_ratio_error_histogram(ratio_error_histogram)

        if args.postFit:
            (ratio_histogram_2_sigma, ratio_error_histogram_2_sigma) = make_ratio_histogram(nBins, xMin, xMax, mc_error_histogram_2_sigma, data_histogram, stack_plot)
            ratio_error_histogram_2_sigma = configure_ratio_error_histogram(ratio_error_histogram_2_sigma)
            ratio_error_histogram_2_sigma.SetFillColor(kYellow)

        if not skip_data:
            g_ratio = make_data_ratio_asymmetric_errors(nBins, xMin, xMax, ggg, data_histogram, mc_error_histogram, stack_plot)
            g_ratio = configure_data_ratio_asymmetric_errors(g_ratio)

            ratio_histogram.SetLineColor(kWhite)
            ratio_histogram.SetMarkerColor(kWhite)
            ratio_histogram.SetFillColor(kWhite)
            ## We're really just drawing the ratio histogram to get the proper dimensions / format
            ## for the ratio histogram box.  The real data ratio histogram is g_ratio.  Unfortunately,
            ## g_ratio is a TGraphAsymmErrors, which doesn't take all the same histo options as a TH1
            ratio_histogram.DrawCopy()
            if args.postFit:
                ratio_error_histogram_2_sigma.DrawCopy('e2same')
            ratio_error_histogram.DrawCopy('e2same')
            ratio_histogram.Draw('sameaxis')
            ## asymmetrical poisson errors for data in ratio plot
            g_ratio.Draw('psame')

        else:
            ratio_histogram.SetLineWidth(cosmetics['g ratio line width'])
            ratio_histogram.DrawCopy()

        l = TLine()
        l.DrawLine(xMin, 1., xMax, 1.)

    if config.get('draw cumulative plot'):
        print 'Making a cumulative plot'
        cumulative_histo = make_cumulative_histo(nBins, xMin, xMax, mc_error_histogram, histogram_dictionary[config['draw cumulative plot']+'_nominal'], distribution)
        cumulative_histo.DrawCopy()
        cumulative_histo.Draw('sameaxis')

    output_dir = os.path.join(out_location, lepton_category+'_'+jet_tag_category)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    plot_name = '%s_%s/%s' % (lepton_category, jet_tag_category, distribution)
    if (config['save png']): canvas.SaveAs(out_location+'/'+plot_name+'.png')
    if (config['save pdf']): canvas.SaveAs(out_location+'/'+plot_name+'.pdf')

    if args.web:
        www_base_directory = plot_helper.get_www_base_directory()
        www_plot_directory = os.path.join(www_base_directory, out_location, lepton_category, jet_tag_category)
        plot_helper.copy_to_www_area(out_location, www_plot_directory, plot_name)

    gPad.Close()
    top_canvas.Close()
    bottom_canvas.Close()
    canvas.Close()
## end draw_stack_plot

def get_group_histogram(distribution, systematic, samples_in_group, lepton_category, jet_tag_category, sample_group=''):
    group_histogram = None
    for sample in samples_in_group:
        try:
            sample_histogram = get_histogram(distribution, systematic, sample, lepton_category, jet_tag_category, sample_group)
            #sample_histogram.Sumw2() ##No longer necessary - done in plot_helper -AWB 30/10/14
            if not group_histogram:
                group_histogram = sample_histogram
            else:
                group_histogram.Add(sample_histogram)
        except Exception, e:
            print 'This is an exception in get_group_histogram'
            print e
            continue

    return group_histogram

## Gets a single histogram.
def get_histogram(distribution, systematic, sample, lepton_category, jet_tag_category, sample_group=''):
    histogram = None

    name = '%s_%s' % (distribution, systematic)
    if systematic == 'nominal':
        name = '%s' % distribution
    if plot_helper.is_matching_data_sample(lepton_categories[lepton_category]['data samples'], sample_group):
        if not 'sideband' in sample and not 'Sideband' in sample:
            name = '%s' % distribution
    if args.postFit:
        name = '%s' % sample
        if args.postFitNorms and not 'err_hist' in sample:
            name = '%s_norm' % sample

    file_name = os.path.join(config['input file location'], lepton_category, '%s_%s_%s_%s.root' % (lepton_category, jet_tag_category, sample, config['input file label']))
    if args.postFit:
        file_name = os.path.join(config['postFit file location'], '%s_%s%s.root' % (lepton_categories[lepton_category]['postFit label'][0], jet_tag_category, lepton_categories[lepton_category]['postFit label'][1]) )
    try:
        root_file = TFile(file_name)
    except ReferenceError:
        raise ReferenceError('Problem finding file %s, skipping it...' % file_name)
    try:
        histogram = root_file.Get(name).Clone()
    except ReferenceError:
        raise ReferenceError('Problem finding %s in file %s, skipping it...' % (name, file_name))

    histogram.SetDirectory(0) ##Decouples histogram from root file

    ## Set the bin statistical error and cross section error
    if not args.postFit:
        sample_info = plot_helper.SampleInformation(sample)
        for bin in range(1, histogram.GetNbinsX()):
            bin_error_squared = math.pow(histogram.GetBinError(bin), 2)
            if not args.postFitNorms:
                bin_error_squared += math.pow(histogram.GetBinContent(bin)*sample_info.x_section_error_ttV/sample_info.x_section, 2)
            histogram.SetBinError(bin, math.sqrt(bin_error_squared))

    if distributions[distribution][0] == '':
        distributions[distribution][0] = histogram.GetXaxis().GetTitle()
    if distributions[distribution][1]:
        histogram = move_underflow_into_hist(histogram)
    if distributions[distribution][2]:
        histogram = move_overflow_into_hist(histogram)

    if sample_group == 'data':
        histogram.SetLineColor(cosmetics['data line color'])
        histogram.SetMarkerColor(cosmetics['data marker color'])
        histogram.SetFillColor(cosmetics['data fill color'])
        histogram.SetLineWidth(cosmetics['data line width'])
        histogram.SetMarkerStyle(cosmetics['data marker style'])
        histogram.SetMarkerSize(cosmetics['data marker size'])
    elif config.has_key('background samples') and (sample_group in background_samples):
        histogram.SetLineColor(background_samples[sample_group]['color'])
        histogram.SetFillColor(background_samples[sample_group]['color'])
        histogram.SetFillStyle(cosmetics['background fill style'])
    elif config.has_key('signal samples') and (sample_group in signal_samples):
        histogram.SetLineColor(signal_samples[sample_group]['color'])
        histogram.SetFillColor(signal_samples[sample_group]['color'])
        if (signal_samples[sample_group]['stack or line'] == 'line'):
            histogram.SetLineWidth(signal_samples[sample_group]['line width'])
            histogram.SetFillStyle(0)
        elif (signal_samples[sample_group]['stack or line'] == 'stack'):
            histogram.SetFillStyle(cosmetics['background fill style'])
        else:
            raise Exception('Error-- must select either line or stack for signal_sample '+sample_group+'.  Quitting.')

    return histogram
## end of get_histogram

def make_legend():
    stack_plot_legend = TLegend(cosmetics['legend x1'], cosmetics['legend y1'], cosmetics['legend x2'], cosmetics['legend y2'])
    stack_plot_legend.SetFillColor( cosmetics['legend fill color'])
    stack_plot_legend.SetBorderSize(cosmetics['legend border size'])
    stack_plot_legend.SetNColumns(cosmetics['legend n columns'])

    return stack_plot_legend
## end of make_legend

def make_info_tex_objects(lepton_category, jet_tag_category):
    luminosity_info_tex = TLatex()
    luminosity_info_tex.SetTextAlign(11)
    luminosity_info_tex.SetNDC()
    luminosity_info_tex.SetTextFont(cosmetics['lumi text font'])
    luminosity_info_tex.SetTextSize(cosmetics['lumi text size'])
    luminosity_info_tex.SetTitle(config['lumi era string'])

    selection_info_tex = TLatex()
    selection_info_tex.SetTextAlign(11)
    selection_info_tex.SetNDC()
    selection_info_tex.SetTextFont(cosmetics['selection text font'])
    selection_info_tex.SetTextSize(cosmetics['selection text size'])
    selection_info_tex.SetTitle('%s%s' % (lepton_categories[lepton_category]['tex name'], jet_tag_categories[jet_tag_category]))
        
    CMS_info_tex = TLatex()
    CMS_info_tex.SetTextAlign(11)
    CMS_info_tex.SetNDC()
    CMS_info_tex.SetTextFont(cosmetics['CMS text font'])
    CMS_info_tex.SetTextSize(cosmetics['CMS text size'])
    CMS_info_tex.SetTitle(config['CMS string'])
    if args.blind:
        CMS_info_tex.SetTitle(config['CMS string'].replace('Preliminary', 'Simulation'))

    fit_status_tex = TLatex()
    fit_status_tex.SetTextAlign(11)
    fit_status_tex.SetNDC()
    fit_status_tex.SetTextFont(cosmetics['fit status text font'])
    fit_status_tex.SetTextSize(cosmetics['fit status text size'])
    fit_status_tex.SetTitle('#it{Pre-fit}')
    if args.postFit:
        fit_status_tex.SetTitle('#it{Post-fit}')
    if args.postFitNorms:
        fit_status_tex.SetTitle('#it{Post-fit yields}')

    SF_info_tex = TLatex()
    SF_info_tex.SetNDC()
    SF_info_tex.SetTextFont(cosmetics['SF text font'])
    SF_info_tex.SetTextSize(cosmetics['SF text size'])
    SF_strings = ['%s x %s' % (signal_samples[sample_group]['draw name'], signal_samples[sample_group]['scale']) for sample_group in signal_samples]
    SF_info_tex.SetTitle(', '.join(SF_strings))

    return luminosity_info_tex, selection_info_tex, CMS_info_tex, fit_status_tex, SF_info_tex
## end of make_info_tex_objects

def move_extra_into_hist(histogram):
    '''Move over/underflow bins into the histogram'''
    num_bins = histogram.GetNbinsX()
    histogram.SetBinContent(1, histogram.GetBinContent(0) + histogram.GetBinContent(1))
    histogram.SetBinContent(0,0)
    histogram.SetBinContent(num_bins, histogram.GetBinContent(num_bins) + histogram.GetBinContent(num_bins+1))
    histogram.SetBinContent(num_bins+1, 0)

    return histogram
## end move_extra_into_histogram()

def move_underflow_into_hist(histogram):
    '''Move underflow bins into the histogram'''
    num_bins = histogram.GetNbinsX()
    histogram.SetBinContent(1, histogram.GetBinContent(0) + histogram.GetBinContent(1))
    histogram.SetBinContent(0,0)

    return histogram
## end move_underflow_into_histogram()

def move_overflow_into_hist(histogram):
    '''Move overflow bins into the histogram'''
    num_bins = histogram.GetNbinsX()
    histogram.SetBinContent(num_bins, histogram.GetBinContent(num_bins) + histogram.GetBinContent(num_bins+1))
    histogram.SetBinContent(num_bins+1, 0)

    return histogram
## end move_overflow_into_histogram()

def configure_stack(stack_plot, plot_max):
    stack_plot.SetTitle(';;Events')
    legend_factor = (cosmetics['legend y2'] - cosmetics['top canvas y1']) / (cosmetics['legend y1'] - cosmetics['top canvas y1'])
    if config['selection info']:
        legend_factor = (cosmetics['legend y2'] - cosmetics['top canvas y1']) / (cosmetics['CMS text y'] - 0.015 - cosmetics['top canvas y1'])
    stack_plot.SetMinimum(cosmetics['stack minimum'])
    stack_plot.SetMaximum( max(cosmetics['stack lowest maximum'], plot_max*legend_factor) )
    if (config['log scale']):
        stack_plot.SetMinimum(cosmetics['stack minimum log scale'])
        stack_plot.SetMaximum( max(cosmetics['stack lowest maximum log scale'], plot_max*cosmetics['stack maximum factor log scale']) )

    return stack_plot
## end configure_stack

def get_configured_canvases():
    top_canvas = TPad('top_canvas', 'up', cosmetics['top canvas x1'], cosmetics['top canvas y1'], cosmetics['top canvas x2'], cosmetics['top canvas y2'])
    bottom_canvas =  TPad ('bottom_canvas', 'down', cosmetics['bottom canvas x1'], cosmetics['bottom canvas y1'], cosmetics['bottom canvas x2'], cosmetics['bottom canvas y2'])

    top_canvas.SetLeftMargin(cosmetics['top canvas left margin'])
    bottom_canvas.SetLeftMargin(cosmetics['bottom canvas left margin'])

    top_canvas.SetRightMargin(cosmetics['top canvas right margin'])
    bottom_canvas.SetRightMargin(cosmetics['bottom canvas right margin'])

    top_canvas.SetTopMargin(cosmetics['top canvas top margin'])
    bottom_canvas.SetBottomMargin(cosmetics['bottom canvas bottom margin'])
    bottom_canvas.SetTickx()

    top_canvas.Modified()
    bottom_canvas.Modified()

    return top_canvas, bottom_canvas
## end get_configured_canvases

def get_configured_data_asymmetric_errors(data_histogram):
    ggg = TGraphAsymmErrors(data_histogram)
    alpha = 1 - 0.6827
    num_bins = ggg.GetN()
    for bin in range(0, num_bins):
        NN = ggg.GetY()[bin]
        if NN == 0:
            LOW = 0
            ggg.SetPoint(bin, ggg.GetX()[bin], 0.001) ## So you can see points with 0 events
        else:
            LOW =  Math.gamma_quantile(alpha/2, NN, 1.)
        UP =  Math.gamma_quantile_c(alpha/2, NN+1, 1)
        ggg.SetPointEYlow(bin, NN-LOW)
        ggg.SetPointEYhigh(bin, UP-NN)

    ggg.SetLineColor(cosmetics['ggg line color'])
    ggg.SetMarkerStyle(cosmetics['ggg marker style'])
    ggg.SetLineWidth(cosmetics['ggg line width'])
    ggg.SetMarkerSize(cosmetics['ggg marker size'])

    return ggg

def make_cumulative_histo(n_bins, x_min, x_max, bg_histo, signal_histo, distribution):
    cumulative_histo = TH1D('cumulative_histo', '', n_bins, x_min, x_max)
    cumulative_histo.SetStats(kFALSE)
    ratio_title_string = ';%s;%s' % (distributions[distribution][0], 'C. BG/sig')
    cumulative_histo.SetTitle(ratio_title_string)
    cumulative_histo.GetYaxis().SetTitleSize(cosmetics['ratio hist y axis title size'])
    cumulative_histo.GetYaxis().SetTitleOffset(cosmetics['ratio hist y axis title offset'])
    cumulative_histo.GetYaxis().CenterTitle()
    cumulative_histo.GetYaxis().SetLabelSize(cosmetics['ratio hist y axis label size'])
    cumulative_histo.GetYaxis().SetNdivisions(cosmetics['ratio hist y axis Ndivisions'])
    cumulative_histo.GetXaxis().SetLabelSize(cosmetics['ratio hist x axis label size'])
    cumulative_histo.GetXaxis().SetLabelOffset(cosmetics['ratio hist x axis label offset'])
    cumulative_histo.GetXaxis().SetTitleOffset(cosmetics['ratio hist x axis title offset'])
    cumulative_histo.GetXaxis().SetTitleSize(cosmetics['ratio hist x axis title size'])

    bg_total = 0
    signal_total = 0
    bg_integral = bg_histo.Integral()
    signal_integral = signal_histo.Integral()
    for i in range(1, n_bins+1):
        bg_total += bg_histo.GetBinContent(i)
        signal_total += signal_histo.GetBinContent(i)
        if (bg_total != 0 and signal_total != 0):
            ratio = (bg_total/bg_integral) / (signal_total/signal_integral)
            cumulative_histo.SetBinContent(i, ratio)

    return cumulative_histo

def make_ratio_histogram(nBins, xMin, xMax, mc_error_histogram, data_histogram, stack_plot):
    #ratio_histogram = TH1D('ratio_histogram', '', nBins, xMin, xMax)
    #ratio_error_histogram = TH1D('ratio_error_histogram', '', nBins, xMin, xMax)

    ratio_histogram = data_histogram.Clone('ratio_histogram')
    ratio_error_histogram = data_histogram.Clone('ratio_error_histogram')
    ratio_histogram.SetTitle('')
    ratio_error_histogram.SetTitle('')

    for i in range(1, nBins+1):
        ## May be different from the value of of mc_error_histogram.GetBinContent(i)
        mc_value = stack_plot.GetStack().Last().GetBinContent(i)
        mc_error = mc_error_histogram.GetBinError(i)
        mc_error_center = mc_error_histogram.GetBinContent(i)
        if args.postFit and abs(1 - mc_error_center/mc_value) > 0.1:
            print 'WARNING: mc_error and stack differ by %0.2f percent in bin %d' % (abs(1 - mc_error_center/mc_value)*100, i)
        data_value = data_histogram.GetBinContent(i)
        data_error = data_histogram.GetBinError(i)
        if (mc_value !=0 and data_value != 0):
            ratio_value = data_value / mc_value
            ratio_error = ratio_value * data_error / data_value
            ratio_error_mc = mc_error / mc_value
            ratio_value = min(ratio_value, 2)
        elif (mc_value == 0 and data_value != 0):
            ratio_value = 2
            ratio_error = 2
        elif (mc_value != 0 and data_value == 0):
            ratio_error_mc = mc_error / mc_value
        if data_value != 0:
            ratio_histogram.SetBinContent(i, ratio_value)
            ratio_histogram.SetBinError(i, ratio_error)
        if mc_value != 0:
            ratio_error_histogram.SetBinContent(i, mc_error_center/mc_value)
            ratio_error_histogram.SetBinError(i, ratio_error_mc)

    return ratio_histogram, ratio_error_histogram
## end make_ratio_histogram

def configure_ratio_histogram(ratio_histogram, distribution):
    ratio_histogram.SetStats(kFALSE)
    ratio_histogram.SetMinimum(cosmetics['ratio hist min'])
    ratio_histogram.SetMaximum(cosmetics['ratio hist max'])
    ratio_title_string = ';%s;%s' % (distributions[distribution][0], 'Data/Pred.')
    ratio_histogram.SetTitle(ratio_title_string)
    
    ratio_histogram.GetYaxis().CenterTitle()
    ratio_histogram.GetYaxis().SetTitleSize(cosmetics['ratio hist y axis title size'])
    ratio_histogram.GetYaxis().SetTitleOffset(cosmetics['ratio hist y axis title offset'])
    ratio_histogram.GetYaxis().SetNdivisions(cosmetics['ratio hist y axis Ndivisions'])
    ratio_histogram.GetYaxis().SetLabelSize(cosmetics['ratio hist y axis label size'])
    ratio_histogram.GetYaxis().SetLabelOffset(cosmetics['ratio hist y axis label offset'])
    
    ratio_histogram.GetXaxis().SetTitleSize(cosmetics['ratio hist x axis title size'])
    ratio_histogram.GetXaxis().SetTitleOffset(cosmetics['ratio hist x axis title offset'])
    ratio_histogram.GetXaxis().SetNdivisions(cosmetics['ratio hist x axis Ndivisions'])
    ratio_histogram.GetXaxis().SetLabelSize(cosmetics['ratio hist x axis label size'])
    ratio_histogram.GetXaxis().SetLabelOffset(cosmetics['ratio hist x axis label offset'])
    
    ratio_histogram.SetLineColor(cosmetics['ratio hist line color'])
    ratio_histogram.SetMarkerColor(cosmetics['ratio hist marker color'])

    if args.blind:
        ratio_histogram.SetMinimum(0)
        ratio_histogram.SetMaximum(0.3)
        ratio_title_string = ';%s;%s' % (distributions[distribution][0], 'Signal/Pred.')
        ratio_histogram.SetTitle(ratio_title_string)
        ratio_histogram.SetLineColor(config['signal samples'][config['signal samples'].keys()[0]]['color'])
        ratio_histogram.SetMarkerColor(config['signal samples'][config['signal samples'].keys()[0]]['color'])

    return ratio_histogram
## end configure_ratio_histogram

def configure_ratio_error_histogram(ratio_error_histogram):
    ratio_error_histogram.SetMarkerColor(cosmetics['ratio err hist marker color']) ## Does nothing
    ratio_error_histogram.SetMarkerSize(0)
    ratio_error_histogram.SetFillColor(cosmetics['ratio err hist fill color'])

    return ratio_error_histogram

def make_data_ratio_asymmetric_errors(nBins, xMin, xMax, ggg, data_histogram, mc_error_histogram, stack_plot):
    ratio_max = cosmetics['ratio hist max']
    g_ratio = TGraphAsymmErrors(ggg.GetN())
    for bin in range(0, g_ratio.GetN()):
        x_point = data_histogram.GetBinCenter(bin+1)
        x_width = 0.5*data_histogram.GetBinWidth(bin+1)
        yG = ggg.GetY()[bin]
        yG_low  = ggg.GetEYlow()[bin]
        yG_high = ggg.GetEYhigh()[bin]
        if not args.blind:
            y_data = data_histogram.GetBinContent(bin+1)
        else:
            y_data = mc_error_histogram.GetBinContent(bin+1)

        ## May be different from the value of of mc_error_histogram.GetBinContent()
        yMC = stack_plot.GetStack().Last().GetBinContent(bin+1)

        if yMC > 0.0:
            yG_ratio = yG/yMC
            yG_ratio_low = yG_low/yMC
            yG_ratio_high = yG_high/yMC
        else:
            yG_ratio = 0.0
            yG_ratio_low = 0.0
            yG_ratio_high = 0.0

        if y_data > -1: ## Include points with 0 data
            g_ratio.SetPoint(bin, x_point, yG_ratio)
            g_ratio.SetPointEYlow(bin, yG_ratio_low)
            g_ratio.SetPointEYhigh(bin, yG_ratio_high)
            g_ratio.SetPointEXlow(bin, x_width)
            g_ratio.SetPointEXhigh(bin, x_width)

            if (yG_ratio > ratio_max and (yG_ratio - yG_ratio_low) < ratio_max):
                minner = yG_ratio_low - (yG_ratio - ratio_max - 0.0001)
                g_ratio.SetPoint(bin, x_point, ratio_max - 0.0001)
                g_ratio.SetPointEYlow(bin, minner)

    return g_ratio
## end make_data_ratio_asymmetric_errors

def configure_data_ratio_asymmetric_errors(g_ratio):
    g_ratio.SetLineColor(cosmetics['g ratio line color'])
    g_ratio.SetLineWidth(cosmetics['g ratio line width'])

    return g_ratio

if __name__ == '__main__':
    main()
