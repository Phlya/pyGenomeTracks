from . GenomeTrack import GenomeTrack
import numpy as np

class LoopTrack(GenomeTrack):
    SUPPORTED_ENDINGS = ['.bedpe']  # this is used by make_tracks_file to guess the type of track based on file name
    TRACK_TYPE = 'loops'
    OPTIONS_TXT = """
	height = 3
	title = Loops
	color = blue
	line width = 1
	"""
    def __init__(self, *args, **kwarg):
        super(LoopTrack, self).__init__(*args, **kwarg)
        if 'line width' not in self.properties:
            self.properties['line width'] = 1
        if 'line style' not in self.properties:
            self.properties['line style'] = 'solid'
        self.max_height = None
        valid_intervals = 0
        intervals = []
        line_number = 0
        with open(self.properties['file'], 'r') as file_h:
            for line in file_h.readlines():
                line_number += 1
                if line.startswith(('browser', 'track', '#', 'chrom1')):
                    continue
                try:
                    chrom1, start1, end1, chrom2, start2, end2 = line.strip().split('\t')[:6]
                except Exception as detail:
                    self.log.error('File not valid. The format is chrom1 start1, end1, '
                                   'chrom2, start2, end2\nError: {}\n in line\n {}'.format(detail, line))
                    exit(1)

                try:
                    start1 = int(start1)
                    end1 = int(end1)
                    start2 = int(start2)
                    end2 = int(end2)
                except ValueError as detail:
                    self.log.error("Error reading line: {}. One of the fields is not "
                                   "an integer.\nError message: {}".format(line_number, detail))
                    exit(1)

                assert start1 <= end1, "Error in line #{}, end1 smaller than start1 in {}".format(line_number, line)
                assert start2 <= end2, "Error in line #{}, end2 smaller than start2 in {}".format(line_number, line)

                if chrom1 != chrom2:
                    self.log.warn("Only loops in same chromosome are used. Skipping line\n{}\n".format(line))
                    continue

                if start2 < start1:
                    start1, start2 = start2, start1
                    end1, end2 = end2, end1

                # each interval spans from the smallest start to the largest end
                intervals.append([chrom1, start1, end1, start2, end2])
                valid_intervals += 1

        if valid_intervals == 0:
            self.log.warn("No valid intervals were found in file {}".format(self.properties['file']))
        intervals = pd.DataFrame(intervals, columns=['chrom', 'start1', 'end1', 'start2', 'end2'])
        file_h.close()
        self.loops = intervals

        if 'color' not in self.properties:
            self.properties['color'] = 'blue'

        if 'alpha' not in self.properties:
            self.properties['alpha'] = 0.8
    
    def plot_y_axis(self, ax, plot_ax):
        pass
    
    def plot_loops(self, ax, loop):
        from matplotlib.patches import Polygon
        width1 = loop[1] - loop[0]
        width2 = loop[3] - loop[2]
        x0 = (loop[1]+loop[2])/2
        y0 = loop[2] - loop[1]
        
        x1 = x0 + width2/2
        y1 = y0 + width2#/2
        
        x2 = (loop[1]+loop[2])/2
        y2 = loop[3] - loop[0]
        
        x3 = x0 - width1/2
        y3 = y0 + width1#/2

        rectangle = Polygon(np.array([[x0, y0], [x1, y1], [x2, y2], [x3, y3]]),
                           facecolor='none', edgecolor=self.properties['color'],
                           linewidth=self.line_width,
                           ls=self.properties['line style'])
        ax.add_artist(rectangle)

    def plot(self, ax, chrom_region, region_start, region_end):
        """
        Makes a rectangle highlighting an interaction between two regions on a
        linear scale representing interactions between Hi-C bins.
        :param ax: matplotlib axis
        :param label_ax: matplotlib axis for labels
        """
        self.max_height = 0
        count = 0
        
        loops_in_region = self.loops.query("(chrom=='%s')&(end1>%s)&(end1<%s)&(start2>%s)&(start2<%s)" %\
                                           (chrom_region, region_start, region_end, region_start, region_end))
        loops_in_region = loops_in_region[['start1', 'end1', 'start2', 'end2']]
        
        for idx, interval in loops_in_region.iterrows():
            if 'line width' in self.properties:
                self.line_width = float(self.properties['line width'])
            else:
                self.line_width = 0.5 * np.sqrt(interval.data)
            
            self.plot_loops(ax, interval)
            count += 1

        self.log.debug("{} loops were plotted".format(count))
        self.log.debug('title is {}'.format(self.properties['title']))
