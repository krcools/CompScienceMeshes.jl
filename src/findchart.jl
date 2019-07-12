function findchart(charts, chart_tree, point)

    for box in boxes(chart_tree, (c,s)->boxesoverlap(c,s,point,0))
        for i in box
            isinclosure(charts[i], point) && return i
        end
    end

    return nothing
end
