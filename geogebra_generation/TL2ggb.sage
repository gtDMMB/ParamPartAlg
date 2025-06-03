from lxml.etree import Element,SubElement,tostring
from sage.geometry.polyhedron.plot import cyclic_sort_vertices_2d
from sage.graphs.graph_coloring import first_coloring
from zipfile import ZipFile, ZIP_DEFLATED
from pathlib import Path
from os.path import abspath, basename, dirname, splitext
import html as html_module
TEMPLATE_PREFIX='ggb_template'

with open('{0}/ggb_prefix.xml'.format(TEMPLATE_PREFIX),'r') as prefix_source:
    xml_prefix = prefix_source.read()

xml_suffix = """</geogebra>"""

# From https://stackoverflow.com/questions/701704/convert-html-entities-to-unicode-and-vice-versa
def unescaped_html(escaped):
    unescaped = html_module.unescape(escaped)
    return unescaped

def rnapoly2xmlstr(P,name):
    xmlCons = xmlRNAPolyConstruction(P,name)
    return "\n".join([xml_prefix,str(xmlCons),xml_suffix])

def rnapoly2ggb(P,name,outdir):
    xml_str = rnapoly2xmlstr(P,name)

    with open('xml_tmp/geogebra.xml','w') as xml_out:
        xml_out.write(xml_str)
    with ZipFile('{0}/{1}.ggb'.format(outdir,name), 'w', ZIP_DEFLATED) as ggb:
        P = Path('xml_tmp/geogebra.xml')
        ggb.write(P,P.name)
        P = Path('{0}/geogebra_defaults2d.xml'.format(TEMPLATE_PREFIX))
        ggb.write(P,P.name)
        P = Path('{0}/geogebra_defaults3d.xml'.format(TEMPLATE_PREFIX))
        ggb.write(P,P.name)
        P = Path('{0}/geogebra_javascript.js'.format(TEMPLATE_PREFIX))
        ggb.write(P,P.name)
        P = Path('{0}/geogebra_thumbnail.png'.format(TEMPLATE_PREFIX))
        ggb.write(P,P.name)

def cyclically_arranged_vxs(polygon):
    return [tuple(v) for v in cyclic_sort_vertices_2d(polygon.vertices())]


# Given a Color it returns a triple [r,g,b] where r,g,b are ints between
# 0 and 255.
Color2rgb = lambda c: [int('0x'+c.html_color()[2*i+1:2*i+3],16) for i in range(3)]

rgb2html = lambda rgb: '#'+''.join([hex(x).zfill(2) for x in rgb])

rgb2Color = lambda rgb: Color(rgb2html(rgb))

def generate_new_label(L,prefix='X'):
    i=1
    while True:
        proposal = '{0}{1:05d}'.format(prefix,i)
        if proposal not in L:
            break
        i=i+1
    return proposal

class xmlObj:
    def __init__(self, *args, **kwargs):
        self.element = Element(*args, **kwargs)

    def __str__(self):
        return tostring(self.element,pretty_print=True,encoding=str)
    def pprint(self):
        print(str(self))

class xmlPoint(xmlObj):
    labels=[]
    def __init__(self,x,y,z=1.0,show='true',label=None):
        self.X=float(x)
        self.Y=float(y)
        self.Z=float(z)
        if label is None:
            self.label = generate_new_label(self.labels,'pt')
            self.labels.append(self.label)
        xmlObj.__init__(self,'element',type='point',label=self.label)

        SubElement(self.element,'show',object=show,label='false')
        SubElement(self.element,'objColor',r='0',g='0',b='0',alpha='0.0')
        SubElement(self.element,'layer',val='0')
        SubElement(self.element,'labelMode',val='0')
        SubElement(self.element, 'animation',step='1',speed='1',type='1',playing='false')
        SubElement(self.element, 'coords', x='{0:.5f}'.format(self.X),
                   y='{0:.5f}'.format(self.Y), z='{0:.5f}'.format(self.Z))
        SubElement(self.element, 'pointSize', val='2')
        SubElement(self.element, 'pointStyle', val='0')

class xmlSegment(xmlObj):
    labels=[]
    def __init__(self, pt1, pt2, label=None):
        if label is None:
            self.label = generate_new_label(self.labels,'S')
            self.labels.append(self.label)
        xmlObj.__init__(self,'element',type='segment',label=self.label)
        (x1,y1) = pt1
        (x2,y2) = pt2
        # Homogeneous coordinate representation for the line through this segment
        # Ax+By+C=0
        self.A = float(y1-y2)
        self.B = float(x2-x1)
        self.C = -self.A*float(x1)-self.B*float(y1)

        SubElement(self.element,'show',object='false',label='false')
        SubElement(self.element,'objColor',r='0',g='0',b='0',alpha='0.0')
        SubElement(self.element,'layer',val='0')
        SubElement(self.element,'labelMode',val='0')
        SubElement(self.element,'auxiliary',val='false')
        SubElement(self.element, 'coords', x='{0:.5f}'.format(self.A),
                   y='{0:.5f}'.format(self.B), z='{0:.5f}'.format(self.C))
        SubElement(self.element, 'lineStyle', thickness='1', type='10', typeHidden='1', opacity='5')
        SubElement(self.element, 'outlyingIntersections', val='false')
        SubElement(self.element, 'keepTypeOnTransform', val='true')

class xmlPolygonCommand(xmlObj):
    def __init__(self,poly_label, pt_labels,seg_labels):
        xmlObj.__init__(self,'command',name='Polygon')
        inpt = SubElement(self.element, 'input')
        inpt.attrib.update([('a{0}'.format(i),p) for i,p in enumerate(pt_labels)])
        outpt=SubElement(self.element,'output',a0=poly_label)
        outpt.attrib.update([('a{0}'.format(i),p) for i,p in enumerate(seg_labels,1)])


class xmlPolygon(xmlObj):
    labels=[]
    def __init__(self,r,g,b,b_val,label=None):
        if label is None:
            self.label = generate_new_label(self.labels,'P')
            self.labels.append(self.label)
        xmlObj.__init__(self,'element',type='polygon',label=self.label)
        SubElement(self.element, 'lineStyle', thickness='3', type='10', typeHidden='1', opacity='114')
        SubElement(self.element,'show',object='true',label='false', ev='4')
        SubElement(self.element, 'condition', showObject= unescaped_html("b &#x225F; {0}".format(int(b_val)) ) )
        SubElement(self.element,'objColor',r='{0}'.format(r),g='{0}'.format(g),b='{0}'.format(b),alpha='0.8')
        SubElement(self.element,'layer',val='0')
        SubElement(self.element,'labelMode',val='0')
        SubElement(self.element,'ggbscript',val='')

class xmlSlider(xmlObj):
    def __init__(self,mini,maxi,label,xpos,ypos):
        xmlObj.__init__(self, 'element', type='numeric', label=label)
        SubElement(self.element, 'value', val='1.0')
        SubElement(self.element, 'show', object='true', label='true')
        SubElement(self.element, 'objColor', r='0',g='0',b='0',alpha='0.1')
        SubElement(self.element, 'layer', val='0')
        SubElement(self.element, 'labelMode', val='1')
        SubElement(self.element, 'slider', min='{0}'.format(int(mini)), max='{0}'.format(int(maxi)),
                  absoluteScreenLocation='true', width='200.0',x='{0}'.format(xpos), y='{0}'.format(ypos),fixed='true',horizontal='true',
                  showAlgebra='true')
        SubElement(self.element, 'lineStyle',thickness='10',type='0',typeHidden='1')
        SubElement(self.element, 'animation',step='1',speed='1',type='0',playing='false')

class xmlExpression(xmlObj):
    def __init__(self,expr,label):
        xmlObj.__init__(self,'expression',label=label, exp=expr)


class xmlText(xmlObj):
    labels=[]
    def __init__(self,x,y,z=1.0,b_val=None,label=None,fixed=False,ltype=None):
        if label is None:
            self.label = generate_new_label(self.labels,'text')
            self.labels.append(self.label)
        xmlObj.__init__(self,'element',type='text',label=self.label)

        self.x = float(x)
        self.y = float(y)
        self.z = float(z)

        SubElement(self.element,'show',object='true',label='true', ev='40')
        if b_val is not None:
            if ltype is None:
                SubElement(self.element,'condition', showObject= unescaped_html('b &#x225F; {0}'.format(int(b_val)) ))
            else:
                SubElement(self.element,'condition', showObject= unescaped_html('b &#x225F; {0} &#x2227; regionLabel &#x225F; {1}'.format(int(b_val),int(ltype)) ))
        SubElement(self.element,'objColor',r='0',g='0',b='0',alpha='0.0')
        SubElement(self.element,'layer',val='0')
        SubElement(self.element,'labelMode',val='0')
        if fixed:
            SubElement(self.element, 'fixed', val='true')
            SubElement(self.element, 'absoluteScreenLocation', x='{0}'.format(int(self.x)), y='{0}'.format(int(self.y)))
        else:
            SubElement(self.element,'startPoint', x='{0:.5f}'.format(self.x),
                   y='{0:.5f}'.format(self.y), z='{0:.5f}'.format(self.z))

class xmlRNAPolyConstruction(xmlObj):
    def __init__(self,rnapoly,poly_name,b_min=0,b_max=0,b_step=10):
        xmlObj.__init__(self,'construction', title=poly_name, author='RNApoly2ggb', date='')

        self.rnapoly = rnapoly

        self.data = colored_layer_graphs(rnapoly,1,b_min,b_max,b_step=0)
        print("Generating points...")
        self.generate_points()
        print("Generating polygons...")
        self.generate_polygons()
        print("Generating labels...")
        self.generate_labels()
        print("Generating sliders...")
        self.bslider = xmlSlider(b_min,b_max,'b',xpos=90.0,ypos=150.0)
        self.lslider = xmlSlider(0,2,'regionLabel', xpos=300.0,ypos=150.0)
        print("Generating global labels...")
        # slices = rnapoly.d1_slices()
        # for a_slice in slices:
        #     break
        # seq_length = len(a_slice.original_structure)
        seq_length = "LENGTH"
        self.meta_label = xmlText(x=90,y=50,fixed=True)
        self.meta_label_expr = xmlExpression('"{0}; len: {1}"'.format(poly_name,seq_length),self.meta_label.label)

        self.append_elements()

    def generate_points(self):
        pt_list=[]
        for b in self.data:
            all_b_regs = self.data[b][0][0]
            for sign,R in all_b_regs:
                #for p in R.vertices_list():
                #    pt_list.append((p[0],p[1],b))
                pt_list.extend([tuple(p) for p in R.vertices_list()])

        self.pt_list = list(set(pt_list))
        self.pt_list.sort()

        self.xmlPoint = { pt:xmlPoint(pt[0],pt[1],show='false') for pt in self.pt_list }

    def generate_polygons(self):

        #self.color_regions()

        self.xmlPolygon = {}
        self.xmlPolygonCommand = {}
        self.region_data = {}
        for b_val in self.data:
            all_b_regs = self.data[b_val][0][0]

            poly_counter = 0
            poly_colors=[]
            for sign,R in all_b_regs:
                poly_counter+=1
                self.region_data[(R,b_val)]={}
                r,g,b = Color2rgb(R.color)
                poly_colors.append((r,g,b))
                xPo = xmlPolygon(r,g,b,b_val)
                self.xmlPolygon[(R,b_val)] = xPo

                vxs = cyclically_arranged_vxs(R)

                regionSegs = [(vxs[i-1],vxs[i]) for i in range(len(vxs))]
                regionXmlSeg = { s:xmlSegment(s[0],s[1]) for s in regionSegs }

                xPolyCo = xmlPolygonCommand(xPo.label, [self.xmlPoint[p].label for p in vxs], [regionXmlSeg[s].label for s in regionXmlSeg] )

                self.xmlPolygonCommand[(R,b_val)]= xPolyCo
                self.region_data[(R,b_val)] = [regionSegs, regionXmlSeg ]

    def color_regions(self):
        # Color the slices using a graph coloring algorithm
        regions = [R for sig,R in self.all]
        G = region_graph(regions)

        # Colorblind-friendly palette from http://bconnelly.net/2013/10/creating-colorblind-friendly-figures/
        plot_colors = (Color(0.9, 0.6, 0), Color(0.35, 0.7, 0.9), Color(0, 0.6, 0.5), Color(0.95, 0.9, 0.25), Color(0, 0.45, 0.7), Color(0.8, 0.4, 0), Color(0.8, 0.6, 0.7))
        # rgb_colors = [[165,0,38],[215,48,39],[244,109,67],[253,174,97],[254,224,144],[255,255,191],[224,243,248],[171,217,233],[116,173,209],[69,117,180],[49,54,149]]
        # plot_colors = [rgb2Color(rgb) for rgb in rgb_colors]
        n_colors = len(plot_colors)
        region_color_list = first_coloring(G)#, n_colors)
        region_color_dict = {region: plot_colors[i] for i in range(len(region_color_list)) for region in region_color_list[i]}

        # Store the colors
        for region in regions:
            region.color = region_color_dict[region]


    def generate_labels(self):
        self.region_label_data = {}
        for b_val in self.data:
            all_b_regs = self.data[b_val][0][0]
            for sig,R in all_b_regs:

                type1text = '"({0},{1})"'.format(sig[0],sig[2])
                type2text = '"({0},{1},{2},{3})"'.format(sig[0],sig[1],sig[2],sig[3])
                x,y = R.center()

                xmlLabel1 = xmlText(x,y,b_val=b_val,ltype=1)
                xmlExpr1 = xmlExpression(type1text,xmlLabel1.label)

                xmlLabel2 = xmlText(x,y,b_val=b_val,ltype=2)
                xmlExpr2 = xmlExpression(type2text,xmlLabel2.label)

                self.region_label_data[(R,b_val)] = ((xmlExpr1,xmlExpr2), (xmlLabel1,xmlLabel2))


    def append_elements(self):
        for pt in self.pt_list:
            self.element.append(self.xmlPoint[pt].element)
        for b_val in self.data:
            all_b_regs = self.data[b_val][0][0]

            for sig,R in all_b_regs:

                self.element.append(self.xmlPolygonCommand[(R,b_val)].element)
                self.element.append(self.xmlPolygon[(R,b_val)].element)
                regionSegs,regionXmlSeg = self.region_data[(R,b_val)]

                for s in regionSegs:
                    self.element.append(regionXmlSeg[s].element)
            for sig,R in all_b_regs:
                xExprs,xLbls = self.region_label_data[(R,b_val)]
                self.element.append(xExprs[0].element)
                self.element.append(xExprs[1].element)
                self.element.append(xLbls[0].element)
                self.element.append(xLbls[1].element)

        self.element.append(self.bslider.element)
        self.element.append(self.lslider.element)

        self.element.append(self.meta_label_expr.element)
        self.element.append(self.meta_label.element)


def cell_graph(P, boundary_dim=1, lower_b = -3, upper_b= 3,details=False):
    import itertools
    from collections import defaultdict

    regions = set()
    slices = []

    for b in range(lower_b,upper_b+1):
        b_regions = sliced_cones_b(P,b)
        slices.append(b_regions)
        regions=regions.union(b_regions)

    edge_dict = {region:[] for region in regions}
    for pair in itertools.combinations(regions, 2):
        a, b = pair
        if a != b and a.intersection(b).dimension()>=boundary_dim:
            edge_dict[a].append(b)

    result = Graph(edge_dict)
    return result if not details else (result,slices)

def classify_regions(bounded_regs):
    all_regs = []
    wedges=[]
    stripes=[]
    polygons=[]
    for br in bounded_regs:
        # s will be the slice from which br was constructed
        s = br.original_region

        # discriminant will be the number of rays defining the slice
        # cones are bounded by 2 rays, stripes by 1 and polygons by 0.
        discr = len(s.rays())

        # We take a point in the corresponding polytope to be the signature
        signature = tuple(br.original_vertex)#cone2poly(s.original_cone).representative_point()

        if discr == 2:
            wedges.append((signature,br))
        elif discr==1:
            stripes.append((signature,br))
        elif discr==0:
            polygons.append((signature,br))
        all_regs.append((signature,br))
    return all_regs,wedges,stripes,polygons

def bounded_regions(regions, points_to_include, xbounds = None, ybounds = None):
    """
    Plot a collection of cone slices in R2
    """
    import itertools
    vertices = [vert for region in regions for vert in region.vertices()] # Find all the vertices to help us construct the bounding box
    vertices.extend(points_to_include) # Append any marked points to ensure they are visible

    # Identify the maximum values of x and y for the box
    xvals = [abs(vert[0]) for vert in vertices]
    xmax = max(xvals)
    xmin = -xmax

    yvals = [abs(vert[1]) for vert in vertices]
    ymax = max(yvals)
    ymin = -ymax

    auto_box_scale = 3/2
    if xbounds is None:
        xbounds = [auto_box_scale * xmin, auto_box_scale * xmax]

    if ybounds is None:
        ybounds = [auto_box_scale * ymin, auto_box_scale * ymax]

    bounding_box = Polyhedron(itertools.product(xbounds, ybounds))

    bounded_regions = set()
    for region in regions:
        bounded_region = bounding_box.intersection(region)

        # Discard the region if it's empty
        if bounded_region.dimension() < 2:
            continue

        # Otherwise, add some metadata and insert it in the result set
        bounded_region.original_region = region
        bounded_region.original_vertex = region.original_vertex
        bounded_region.original_structure = region.original_structure
        bounded_regions.add(bounded_region)
        
    return bounded_regions

#@time_wrap
def colored_layer_graphs(TLs, boundary_dim=1, lower_b=0, upper_b=0.5, b_step=0.1, x=0,X=2048,y=-2048,Y=2048,debug=False):
    import itertools
    # Colorblind-friendly palette from http://bconnelly.net/2013/10/creating-colorblind-friendly-figures/
    plot_colors = (Color(0.9, 0.6, 0), Color(0.35, 0.7, 0.9), Color(0, 0.6, 0.5), Color(0.95, 0.9, 0.25), Color(0, 0.45, 0.7), Color(0.8, 0.4, 0), Color(0.8, 0.6, 0.7))

    if debug:
        first_coloring = time_wrap(graph_coloring.first_coloring)
    else:
        first_coloring = graph_coloring.first_coloring
    data={b:[] for b in range(lower_b,upper_b+1)}
    graphs=[]
    # For each value of b we compute the corresponding graph
    for b in range(lower_b,upper_b+1):
        regions = TLs[b]
        #regions=regions.union(b_regions)
        bounded_regs = bounded_regions(regions, (), (x,X),(y,Y))
        all_regs,wedges,stripes,polygons = classify_regions(bounded_regs)
        edge_dict = {region.original_vertex:[] for region in bounded_regs}
        for pair in itertools.combinations(bounded_regs, 2):
            u, v = pair
            if u != v and u.intersection(v).dimension()>=boundary_dim:
                edge_dict[u.original_vertex].append(v.original_vertex)

        graphs.append(Graph(edge_dict))
        data[b].append([all_regs,wedges,stripes,polygons])
        data[b].append(graphs[-1])

    # Now we compute the vertices that appear in at least two layers
    all_vxs = set([])
    for g in graphs:
        all_vxs=all_vxs.union(g.vertices())
    # common_vxs = [ v for v in all_vxs if any([v in graphs[i] and v in graphs[i+1] for i in range(len(graphs)-1)])]
    common_vxs = all_vxs
    # We now take the union of the graphs induced by the common vertices to assign them a fixed color class
    # when coloring each of the layers
    #common_edges = list(set([(u,v) for u,v in itertools.combinations(common_vxs, 2) if any([g.has_edge(u,v) for g in graphs])]))
    #common_graph = Graph(common_edges).to_simple()
    common_graph=Graph()
    for g in graphs:
        common_graph=common_graph.union(g.subgraph(common_vxs))

    reg2int = common_graph.relabel(return_map=True)
    int2reg = {reg2int[r]:r for r in reg2int}

    color_classes = first_coloring(common_graph,min(6, len(g.vertices())))
    n_common_colors = len(color_classes)

    temp_coloring = []

    # We iterate through each graph g and color a supergraph h
    # of g. We store the coloring and some extra info to
    # later adjust the coloring using this info.
    for ig,g in enumerate(graphs):
        h=g
        # h is a supergraph of g that contains extra vertices: 1,...,k.
        # These extra vertices form a clique and vx i is adjacent to the
        # neighbors of vertices in color_classes[i]. This is to make sure
        # we can color all vertices in color_classes[i] with
        # whatever color i was assigned in the coloring of h.

        # We begin by adding the vertices and joining it to the neighbors
        # of vxs in color_classes[i].
        for i in range(n_common_colors):
            h.add_vertex(i)
            for v in color_classes[i]:
                if int2reg[v] in h:
                    for u in h.neighbors(int2reg[v]):
                        h.add_edge(i,u)

        # Now we form the clique from the new vertices.
        for i,j in itertools.combinations(range(n_common_colors),2):
            h.add_edge(i,j)

        # We store the coloring of g together with the color that
        # vertex i was assigned, for i=1,...,k
        g_color_classes = first_coloring(h,7)
        i_color ={i:min([j for j,c in enumerate(g_color_classes) if i in c]) for i in range(n_common_colors)}

        temp_coloring.append((g_color_classes,i_color))

    final_coloring = []
    for ig,g in enumerate(graphs):
        b = lower_b+ig
        # Let us now adjust the color classes.
        # For each coloring of g, we will make color assigned
        # to vx i the ith_color. These are the relevant color classes
        # that need to be preserved throught copies of g. All other
        # color classes are only locally relevant.
        new_color_classes = []
        i_color = temp_coloring[ig][1]
        old_color_classes = temp_coloring[ig][0]

        for i in range(n_common_colors):
            # Begin by copying the color class that contains vx i

            ith_color_class = copy(old_color_classes[i_color[i]])
            # we now delete the dummy vx i
            ith_color_class.remove(i)

            n_del = 0
            n_add = 0
            # We now add all vxs in color_classes[i] and delete all vxs
            # in color_classes[j] with i!=j.
            for j in range(n_common_colors):
                if i == j:
                    # Since i is adjacent to all neighbours of vxs in color_classes[i]
                    # then we can just add all of color_classes[i] to the ith color class.
                    # (aka we can rest asured that by adding all vxs in color_classes[i]
                    # to this list there will be no edges within this set)
                    for k in color_classes[i]:
                        v = int2reg[k]
                        if v in g and v not in ith_color_class:
                            ith_color_class.append(v)
                            n_add += 1
                else:
                    # here we delete any common vertex that appears in this
                    # color class and is not in color_classes[i] (these will
                    # reappear in their corresponding ith_color_class).
                    for k in color_classes[j]:
                        v = int2reg[k]
                        if v in g and v in ith_color_class:
                            ith_color_class.remove(v)
                            n_del+=1

            new_color_classes.append(ith_color_class)

        # we are only left to handle the remaining color classes
        # order for these classes doesn't really matter since
        # all of their vertices will only appear in the current
        # graph g
        used_colors = set(i_color.values())
        remaining_colors = set(range(len(old_color_classes))).difference(i_color.values())

        for i in remaining_colors:
            ith_color_class = copy(old_color_classes[i])
            # however, we should make sure no vertex from the common
            # vertices appears in these color classes, so we remove
            # them.
            for v in common_vxs:
                if v in ith_color_class:
                    ith_color_class.remove(v)
            new_color_classes.append(ith_color_class)
        final_coloring.append(new_color_classes)
        data[b].append(new_color_classes)

    for b in data:
        final_color_classes = data[b][2]
        region_color_dict = {region:plot_colors[i] for i in range(len(final_color_classes)) for region in final_color_classes[i]}
        all_regs = data[b][0][0]
        for _,reg in all_regs:
            R=reg.original_vertex
            reg.color = region_color_dict[R]
    
    data[1] = data[0]
    return data

def generate_ggb(infname, db=True):
    outdir = dirname(abspath(infname))		

    if db:
        with open('acc_fname.pkl', 'rb') as pk_file:
            acc_dens_fnames_pairs=pickle.load(pk_file)

        filenames = iterate_sorted('sobjs/','*.sobj')

        for prefix, fname in acc_dens_fnames_pairs:

            code = get_code(fname)
            try:
                filename = get_filename_with_code(filenames,code)
            except IndexError:
                print("[{0}] No sobj file found for {1}.".format(datetime.now(),code))
                print("*"*30)

                continue

            if filename.find(infname) < 0:
                continue

            base_name = splitext(basename(filename))[0]

            print("[{0}] Sarting process for  {1}.".format(datetime.now(),base_name ))

            P = RNAPolytope.construct_from_file(filename)

            try:
                rnapoly2ggb(P,'{0}--{1}'.format(prefix,base_name),outdir)
            except TypeError:
                print("Type error in colored_layer_graphs")
            print("[{0}] Done  {1}.".format(datetime.now(),base_name ))
    else:
        base_name = splitext(basename(infname))[0]
        TL_regions = construct_TL_regions_from_file(infname)
        print(TL_regions)
        rnapoly2ggb(TL_regions,base_name,outdir)

def generate_ggbs():#xbds=(-200,200), ybds=(-200,200)):
    with open('acc_fname.pkl', 'rb') as pk_file:
        acc_dens_fnames_pairs=pickle.load(pk_file)

    filenames = iterate_sorted('sobjs/','*.sobj')

    for prefix, fname in acc_dens_fnames_pairs:

        code = get_code(fname)
        try:
            filename = get_filename_with_code(filenames,code)
        except IndexError:
            print("[{0}] No sobj file found for {1}.".format(datetime.now(),code))
            print("*"*30)

            continue

        base_name = splitext(basename(filename))[0]

        print("[{0}] Sarting process for  {1}.".format(datetime.now(),base_name ))

        P = RNAPolytope.construct_from_file(filename)

        try:
            rnapoly2ggb(P,'{0}--{1}'.format(prefix,base_name))
        except TypeError:
            print("Type error in colored_layer_graphs")
        print("[{0}] Done  {1}.".format(datetime.now(),base_name ))


def construct_TL_regions_from_file(filename, b_val = 0):
    regions = []
    with open(filename) as f:
        for l in f.readlines():
            if len(l) <= 1:
                continue
            
            sig, points = l.strip("\n").split(":")
            sig = eval(sig)
            if len(points) == 0 or "inf" in points:
                continue

            points = tuple(eval(p) for p in points.split(";"))
            r_points = tuple((Rational(a), Rational(c)) for a, c in points)
            region = Polyhedron(r_points)
            region.original_vertex = sig
            region.original_structure = ":)"
            regions.append(region)
    
    b_val_map = {b_val: regions}
    return b_val_map

        
