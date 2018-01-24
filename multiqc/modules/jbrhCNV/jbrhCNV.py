from multiqc import config
from multiqc.plots import scatter,table,linegraph
from multiqc.modules.base_module import BaseMultiqcModule
from collections import OrderedDict
import logging
import re

log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='jbrhCNV',anchor='jbrhCNV',
            href="http://www.jabrehoo.com",info = 'is an tool estimate human '\
            'chromosome copy number varitions based on next generation sequencing data.')

        # set up data structure
        self.cnv_data = {
            'stat':{},
            'report':{},
            'winplot':{},
        }

        # Find and parse cnv stat data
        for f in self.find_log_files('jbrhCNV/stat'):
            parsed_data = self.parse_cnv_stat(f['f'])
            if parsed_data is not None:
                if f['s_name'] in self.cnv_data['stat']:
                    log.debug("Duplicate stat sample log found! Overwriting: {}".format(f['s_name']))
                self.add_data_source(f,section='stat')
                self.cnv_data['stat'][f['s_name']] = parsed_data

        # Find and parse cnv report data
        for f in self.find_log_files('jbrhCNV/report'):
            #sd = self.cnv_data['stat'][f['s_name']]['*SD']
            parsed_data = self.parse_cnv_report(f['f'])
            if parsed_data is not None:
                if f['s_name'] in self.cnv_data['report']:
                    log.debug("Duplicate report sample log found! Overwriting: {}".format(f['s_name']))
                self.add_data_source(f,section='report')
                self.cnv_data['report'][f['s_name']] = parsed_data

        num_parsed = len(self.cnv_data['stat'])
        num_parsed += len(self.cnv_data['report'])

        if num_parsed == 0:
            raise UserWarning

        # Basic Stats table
        self.cnv_stats_table()

        # Write out to the Report
        if len(self.cnv_data['report']) > 0:
            log.info("Found {} jbrhCNV CNV reports".format(len(self.cnv_data['report'])))
            self.cnv_report_table()



        # Find and parse cnv cnt data and plot one by one
        for f in self.find_log_files('jbrhCNV/winplot'):
            parsed_data = self.parse_cnv_winplot(f['f'])
            if parsed_data is not None:
                if f['s_name'] in self.cnv_data['winplot']:
                    log.debug("Duplicate winplot sample log found! Overwriting: {}".format(f['s_name']))
                self.add_data_source(f,section='winplot')
                self.cnv_data['winplot'][f['s_name']] = parsed_data
                self.cnv_winplot_plot(parsed_data,f['s_name'])

        if len(self.cnv_data['winplot']) > 0:
            log.info("Found {} jbrhCNV winplot reports".format(len(self.cnv_data['winplot'])))


    def parse_cnv_stat(self, stat_content):
        """ get total_reads, MT_ratio, map_ratio, Dup, GC and SD """
        parsed_data = {}

        tmp = stat_content.splitlines()
        tmp1 = tmp[1].split()
        tmp2 = tmp[2].split()
        for i in range(len(tmp1)):
            parsed_data[tmp1[i]] = tmp2[i]

        return parsed_data

    def parse_cnv_report(self,report_content):
        """ get report from report.txt"""
        parsed_data = {}
        parsed_data['report'] = report_content

        return parsed_data

    def parse_cnv_winplot(self,cnt):
        """get each chromosome window rcids and plot"""
        data = list()
        n = 0
        rcid_dic = dict()
        rcid_dic_sex = dict()
        color1 = '#0000FF'
        color2 = '#FF0000'
        color = ""
        for line in cnt.splitlines():
            tmp = line.split()
            chro = re.search("chr(\S+)",tmp[0])
            x = re.match(r'[X|Y]',chro.group(1),)
            if x:
                rcid_dic_sex[chro.group(1),int(tmp[1]),int(tmp[2])] = tmp[3]
            else:
                rcid_dic[int(chro.group(1)),int(tmp[1]),int(tmp[2])] = tmp[3]
        for c,s in sorted(rcid_dic.items()):
            n+=1
            name = 'chr' + str(c[0]) + "_" + str(c[1]) + "_" + str(c[2])
            if int(c[0]) % 2 == 0:
                color = color1
            else:
                color = color2
            data.append(
                {'x':n,'y':float(s),'color':color,'name':name}
            )
        for c,s in sorted(rcid_dic_sex.items()):
            n+=1
            name = 'chr' + str(c[0]) + "_" + str(c[1]) + "_" + str(c[2])
            if c[0] == "X":
                color = color2
            if c[0] == "Y":
                color = color1
            data.append(
                {'x':n,'y':float(s),'color':color,'name':name}
            )
        return data


    def cnv_winplot_plot(self,cnt,s_name):
        data = dict()
        config={
            'title':s_name,
            'ymax':40,
            'ymin':0,
            'marker_size':2,
            'marker_line_width':0
        }
        data[s_name] = cnt
        self.add_section(
            name = s_name,
            anchor = 'wp' + s_name,
            content = scatter.plot(data,config),
        )


    def cnv_stats_table(self):
        """ take the parsed stats from the cnv stat and add them to the basic
        stats table at the top of the report """

        headers = {
            'stat':OrderedDict(),
        }
        headers['stat']['Total_Reads'] = {
            'title': 'Total_Reads',
            'description':'total reads',
            'scale':'Greens',
            'format':'{:,.0f}',
            #'modify': lambda x: x / 1000
        }
        headers['stat']['MT_ratio(%)'] = {
            'title': '% MT',
            'description':'mt reads ratio',
            'scale':'Greens',
            'format':'{:,.4f}',
            'max':0.5,
            'min':0
        }
        headers['stat']['Map_Ratio(%)'] = {
            'title': '% Map_Reads',
            'description':'mapping genome reads / total reads',
            'scale':'RdYlGn',
            'max':100,
            'min':50,
        }
        headers['stat']['Duplicate(%)'] = {
            'title': '% Duplicate',
            'description':'duplicate reads / total reads ',
            'scale':'RdYlGn-rev',
            'max':30,
            'min':0,
        }
        headers['stat']['GC_Count(%)'] = {
            'title': '% GC',
            'description':'GC percent in total reads',
            'scale':'RdYlGn-rev',
            'max':60,
            'min':30,
        }
        headers['stat']['*SD'] = {
            'title': 'SD',
            'description':'average standard deviation of chromosome window RCids',
            'scale':'RdYlGn-rev',
            'max':5,
            'min':0,
        }
        self.general_stats_addcols(self.cnv_data['stat'],headers['stat'])

    def cnv_report_table(self):
        headers = {
            'report':OrderedDict(),
        }
        headers['report']['report'] = {
            'title':'Report',
            'description':'chromosome duplicates and deletions report',
        }
        headers['report']['sd'] = {
            'title': 'SD',
            'description':'average standard deviation of chromosome window RCids',
            'scale':'RdYlGn-rev',
            'max':5,
            'min':0,
        }
        self.add_section(
            name = 'report',
            anchor = 'cnv-report',
            content = table.plot(self.cnv_data['report'],headers['report'])
        )
