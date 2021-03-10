from ... import utils
from ... import ops

plt = utils.py.loadExternalModule('matplotlib', 'pylab')
np  = utils.py.loadExternalModule('numpy')

#########################################################

def forest_group(group_name, test_names, tests, conditions, ax,
           col_effect, col_std, col_qvalue, col_test, col_condition,
           group_header_height=0.17, test_height=0.15, condition_marker_size=1, condition_linewidth=1,
           test_colors=['#EEEEEE','white'], condition_colors=['#a6cee3','#1f78b4','#b2df8a','#33a02c'],
           condition_markers=['s', 'd', 's', 'd'],
           x_offset=0, x_end=2, y_offset=0, se_multiplier=1.96):
    """
    parameters:
    -----------
    group_name: String
    test_names: List[String]
    tests: Dataframe
    conditions: list
    col_effect: String
    col_std: String
    col_qvalue: String
    col_test: String
    col_condition: String
    group_header_height: float
    test_height: float
    x_offset: float,
    y_offset: float
    """
    
    from matplotlib.pylab import Rectangle
    
    rect = Rectangle((x_offset,y_offset), width=(x_end-x_offset), height=-group_header_height,
                     edgecolor="none", facecolor="white", linewidth=0, zorder=2)
    ax.add_patch(rect)
    ax.text(x_offset+0.05, y_offset-group_header_height, group_name.upper(),
            horizontalalignment='left', verticalalignment='bottom',
           fontsize='large')
    
    for i, test_name in enumerate(test_names):
        rect_xlim = [ x_offset, x_end ]
        rect_ylim = [ y_offset-group_header_height-i*test_height,
                      y_offset-group_header_height-(i+1)*test_height ]

        rect = Rectangle([rect_xlim[0], rect_ylim[0]], width=(x_end-x_offset), height=-test_height,
                         edgecolor="none", facecolor=test_colors[i%len(test_colors)])
        ax.add_patch(rect)
        ax.text(x_offset+0.1, np.mean(rect_ylim), test_name,
                horizontalalignment='left', verticalalignment='center', fontsize='small')
        
        condition_y_offset_step = test_height / (len(conditions)+1)
        T = tests[(tests[col_test]==test_name)].set_index(col_condition)
        effects = T[col_effect].to_dict()
        stds    = T[col_std].to_dict()
        qvals   = T[col_qvalue].to_dict()
        for j, cond in enumerate(conditions):
            if cond not in effects:
                continue
            #fi
            
            cond_ypos = rect_ylim[0]-(j+1)*condition_y_offset_step
            
            ax.plot([effects[cond]-se_multiplier*stds[cond], effects[cond]+se_multiplier*stds[cond]],
                    [cond_ypos, cond_ypos],
                    c='k', zorder=1, linewidth=condition_linewidth)
            
            ax.scatter([effects[cond]], cond_ypos, zorder=2,
                       c=condition_colors[j%len(condition_colors)],
                       marker=condition_markers[j%len(condition_markers)],
                       s=condition_marker_size, edgecolors='k', linewidth=condition_linewidth)
            
            ax.text(x_end-0.1, cond_ypos, 'p<0.001' if qvals[cond] < 0.001 else '%0.3f' % qvals[cond],
                    horizontalalignment='right', verticalalignment='center', fontsize='xx-small',
                    color=condition_colors[j%len(condition_colors)])
            if qvals[cond] < 0.05:
                ax.text(x_end-0.1, cond_ypos, '*',
                    horizontalalignment='left', verticalalignment='center', fontsize='xx-small',
                    color=condition_colors[j%len(condition_colors)])
        #efor
            
    #efor
    
        

    
#edef

def forest(test_groups, tests, conditions=None,
           col_effect='effect', col_std='std.error', col_qvalue='qvalue', col_test='gene', col_condition='contr',
           group_header_height=0.03, test_height=None,
           test_colors=['#EEEEEE','white'], condition_colors=['#D31D1F','#2C7BB6','#33a02c'],
           condition_marker_size=8, condition_linewidth=0.75,
           condition_markers=['s', 's', 's', 'D', 'D', 'D'],
           test_left_offset=None,
           test_right_offset=1,
           ticks=[-1,-0.5,0,0.5,1], ax=None, se_multiplier=1.96,
           legend_loc='lower left', legend=True):

    
    if conditions is None:
        conditions = sorted(tests[col_condition].unique())
    #fi
    
    tests = tests[tests[col_test].isin(ops.lst.flatten(test_groups.values()))]
    
    if test_left_offset is None:
        test_left_offset = max(2, 0.09 * max([len(str(t)) for t in ops.lst.flatten(test_groups.values())]))
    #fi
    
    if test_height is None:
        test_height = 0.025/2 * len(conditions)
    #fi
    
    forest_xlim = [ min(tests[col_effect] - tests[col_std]) - test_left_offset,
                    max(tests[col_effect] + tests[col_std]) + test_right_offset ]
    

    
    if ax is None:
        x_size = 2+test_left_offset
        y_size = 0.5 + \
                 0.3*len(test_groups) + \
                 10*test_height*len((ops.lst.flatten(test_groups.values())))
                                             
        print(x_size, y_size)
        fig, axes = utils.figure.subplots(ncols=1, figsize=(x_size, y_size))
        ax = axes[0]
    #fi
    
    group_heights = { g : group_header_height + test_height*len(t) for (g,t) in test_groups.items() }

    top      = sum(group_heights.values())
    y_offset = top
    
    for group_name, test_names in test_groups.items():
        forest_group(group_name, test_names, tests, conditions, ax,
               col_effect, col_std, col_qvalue, col_test, col_condition,
               group_header_height=group_header_height, test_height=test_height,
               test_colors=test_colors, condition_colors=condition_colors,
               condition_markers=condition_markers, condition_marker_size=condition_marker_size,
               condition_linewidth=condition_linewidth,
               x_offset=forest_xlim[0], x_end=forest_xlim[1], y_offset=y_offset,
               se_multiplier=se_multiplier)
        y_offset = y_offset - group_heights[group_name]
    #efor
    
    def strip_trailing_zeros(s):
        while s[-1] in ['0', '.']:
            s = s[:-1]
        #ewhile
        return s
    #edef
        
    
    for tick in ticks:
        ax.plot([tick, tick],
                [-test_height/3-0.01, top],
                'k:' if tick != 0 else 'k-', zorder=1, linewidth=0.75)
        ax.text(tick, -test_height/3-0.015, strip_trailing_zeros('%0.2f' % tick) if int(tick) != tick else '%d' % tick,
                horizontalalignment='center', verticalalignment='top', fontsize='x-small')
    #efor

    
    ax.plot([min(ticks)*1.1, max(ticks)*1.1], [-test_height/3, -test_height/3], 'k-', linewidth=0.75)
    
    ax.set_xlim(forest_xlim)
    ax.set_ylim([-test_height/3-0.01, top])
    ax.axis('off')
    
    if legend:
        from matplotlib.lines import Line2D
        legend_elements = [
            Line2D([0], [0],
                   marker=condition_markers[i%len(condition_markers)],
                   color='k',
                   markerfacecolor=condition_colors[i%len(condition_colors)],
                   label=conditions[i], markersize=condition_marker_size)
            for i in range(len(conditions))
        ]
        ax.legend(handles=legend_elements, loc=legend_loc)
    #fi
    
    return ax.get_figure()
#edef