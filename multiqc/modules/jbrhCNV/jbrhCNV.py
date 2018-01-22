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
            sd = self.cnv_data['stat'][f['s_name']]['*SD']
            parsed_data = self.parse_cnv_report(f['f'],sd)
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


    def parse_cnv_stat(self, stat_content):
        """ get total_reads, MT_ratio, map_ratio, Dup, GC and SD """
        parsed_data = {}

        tmp = stat_content.splitlines()
        tmp1 = tmp[1].split()
        tmp2 = tmp[2].split()
        for i in range(len(tmp1)):
            parsed_data[tmp1[i]] = tmp2[i]

        return parsed_data

    def parse_cnv_report(self,report_content,sd):
        """ get report from report.txt"""
        parsed_data = {}

        parsed_data['report'] = report_content
        parsed_data['sd'] = sd

        return parsed_data

    def parse_cnv_winplot(self,cnt):
        """get each chromosome window rcids and plot"""
        data = list()
        n = 0
        for line in cnt.splitlines():
            tmp = line.split()
            n += 1
            chro = re.search("chr(\S+)",tmp[0])
            color1 = '#0000FF'
            color2 = '#FF0000'
            color = ""
            name = tmp[0] + "_" + tmp[1] + "_" + tmp[2]
            if chro.group(1) == "X":
                color = color2
            elif chro.group(1) == "Y":
                color = color1
            else:
                if int(chro.group(1)) % 2 == 0:
                    color = color1
                else:
                    color = color2
            data.append(
                {'x':n,'y':float(tmp[3]),'color':color,'name':name}
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
        data[s_name]=cnt
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
