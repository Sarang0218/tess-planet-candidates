"""Read the official pipeline's own diagnostics from a SPOC / TESS-SPOC data-validation XML report: centroid offsets of the difference image,
odd/even statistic, bootstrap false-alarm, ghost diagnostic, secondary, fit.   Usage: deep2_dvr.py file.xml [period_to_match]"""
import re, sys
s = open(sys.argv[1]).read(); want = float(sys.argv[2]) if len(sys.argv) > 2 else None
def attrs(tag, txt):
    m = re.search(r'<dv:' + tag + r'\b([^>]*)', txt); return dict(re.findall(r'(\w+)="([^"]*)"', m.group(1))) if m else {}
for pr in re.findall(r'<dv:planetResults\b.*?</dv:planetResults>', s, re.S):
    hd = dict(re.findall(r'(\w+)="([^"]*)"', re.search(r'<dv:planetResults\b([^>]*)', pr).group(1)))
    par = {m[0]: (float(m[2]), float(m[1])) for m in re.findall(r'<dv:modelParameter\b[^>]*name="(\w+)"[^>]*uncertainty="([^"]*)"[^>]*value="([^"]*)"', re.search(r'<dv:allTransitsFit\b.*?</dv:allTransitsFit>', pr, re.S).group(0))}
    P = par.get('orbitalPeriodDays', (float('nan'),))[0]
    if want and abs(P / want - 1) > 0.01: continue
    pc = attrs('planetCandidate', pr); fit = attrs('allTransitsFit', pr)
    print(f"planet {hd.get('planetNumber')}: P = {P:.5f} d, depth {par.get('transitDepthPpm', (0, 0))[0]:.0f} ± {par.get('transitDepthPpm', (0, 0))[1]:.0f} ppm, Rp {par.get('planetRadiusEarthRadii', (0, 0))[0]:.2f} ± {par.get('planetRadiusEarthRadii', (0, 0))[1]:.2f} Re, b {par.get('minImpactParameter', (0, 0))[0]:.2f}, model fit SNR {fit.get('modelFitSnr')}")
    print('  detection:', {k: pc.get(k) for k in ('maxMultipleEventSigma', 'maxSingleEventSigma', 'observedTransitCount', 'expectedTransitCount', 'suspectedEclipsingBinary', 'weakSecondaryStruct') if k in pc})
    bs = attrs('bootstrapResults', pr); print('  bootstrap false-alarm:', {k: bs.get(k) for k in bs if 'alse' in k or 'significance' in k.lower()})
    for tag in ('oddEvenTransitDepthComparisonStatistic', 'shorterPeriodComparisonStatistic', 'longerPeriodComparisonStatistic'):
        a = attrs(tag, pr); print(f'  {tag}:', a)
    ws = attrs('weakSecondary', pr); print('  weak secondary:', {k: ws.get(k) for k in ('maxMes', 'maxMesPhaseInDays', 'mesMad') if k in ws})
    for tag in ('coreApertureCorrelationStatistic', 'haloApertureCorrelationStatistic'):
        print(f'  ghost {tag}:', attrs(tag, pr))
    cr = re.search(r'<dv:differenceImageMotionResults\b.*?</dv:differenceImageMotionResults>', pr, re.S)
    if cr:
        for tag in ('msTicCentroidOffsets', 'msControlCentroidOffsets'):
            blk = re.search(r'<dv:' + tag + r'\b.*?</dv:' + tag + '>', cr.group(0), re.S)
            if blk:
                o = {t: attrs(t, blk.group(0)) for t in ('meanRaOffset', 'meanDecOffset', 'meanSkyOffset')}
                print(f"  {tag} (arcsec): RA {float(o['meanRaOffset'].get('value', 'nan')):+.2f} ± {float(o['meanRaOffset'].get('uncertainty', 'nan')):.2f}, Dec {float(o['meanDecOffset'].get('value', 'nan')):+.2f} ± {float(o['meanDecOffset'].get('uncertainty', 'nan')):.2f}, sky {float(o['meanSkyOffset'].get('value', 'nan')):.2f} ± {float(o['meanSkyOffset'].get('uncertainty', 'nan')):.2f}  ({float(o['meanSkyOffset'].get('value', 'nan'))/max(float(o['meanSkyOffset'].get('uncertainty', 'nan')), 1e-9):.1f} sigma)")
    sd = par.get('starDensitySolarDensity'); print('  other fit parameters:', {k: round(v[0], 4) for k, v in par.items() if k in ('ratioSemiMajorAxisToStarRadius', 'transitDurationHours', 'transitIngressTimeHours', 'starDensitySolarDensity', 'equilibriumTempKelvin', 'effectiveStellarFlux')})
