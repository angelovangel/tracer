// Shiny-compatible chromatogram visualization using D3
(function() {
  const BASE_COLORS = {
    A: '#00D100', // Green
    C: '#0000FF', // Blue
    G: '#000000', // Black
    T: '#FF0000'  // Red
  };

  const QV_COLORS = {
    fail: '#F44336',    // Red
    suspect: '#FB7B26', // Amber
    pass: '#808080'     // Green
  };

  // Main render function that Shiny will call
  function renderChromatogram(el, data) {
    // Clear previous content
    d3.select(el).selectAll('*').remove();

    const { traces, bases, peakLocations, qualityScores } = data;
    if (!traces || !Object.keys(traces).length) return;

    // Setup dimensions
    const margin = { top: 40, right: 20, bottom: 30, left: 50 };
    const width = data.plotWidth || el.offsetWidth; // Use specified width or container width
    const height = 250;
    const plotWidth = width - margin.left - margin.right;
    const plotHeight = height - margin.top - margin.bottom;

    // Create SVG
    const svg = d3.select(el)
      .append('svg')
      .attr('width', width)
      .attr('height', height)
      .append('g')
      .attr('transform', `translate(${margin.left},${margin.top})`);

    // Process trace data - one line per channel
    const traceKeys = Object.keys(traces).filter(k => k.startsWith('DATA'));
    const traceData = traceKeys.map(key => {
      const values = traces[key];
      const points = Array.from(values, (y, x) => ({ x, y }));
      return { key, points };
    });

    // Option 2: Shift so first base is at left edge
    let xShift = 0;
    if (peakLocations && peakLocations.length > 0) {
      xShift = peakLocations[0];
    }

    // Setup scales (absolute, no shift)
    const xExtent = [0, Math.max(...traceData[0].points.map(d => d.x))];
    const yExtent = [
      0,
      Math.max(...traceData.map(trace => 
        Math.max(...trace.points.map(d => d.y))
      ))
    ];

    const xScale = d3.scaleLinear()
      .domain(xExtent)
      .range([0, plotWidth]);

    const yScale = d3.scaleLinear()
      .domain(yExtent)
      .range([plotHeight, 0]);

    // Add traces (one per channel)
    const line = d3.line()
      .x(d => xScale(d.x))
      .y(d => yScale(d.y));

    // Map DATA channels to colors (standard ABI mapping)
    const channelToBase = {
      DATA1: 'G', // Black
      DATA2: 'A', // Green
      DATA3: 'T', // Red
      DATA4: 'C'  // Blue
    };

    // Set lower opacity for raw signal view (when no bases provided)
    const isRawSignal = !bases || !peakLocations;
    const strokeOpacity = isRawSignal ? 1.0 : 0.5;
    const strokeWidth = isRawSignal ? 1.5 : 1.5;
    
    traceData.forEach(trace => {
      svg.append('path')
        .datum(trace.points)
        .attr('fill', 'none')
        .attr('stroke', BASE_COLORS[channelToBase[trace.key]])
        .attr('stroke-width', strokeWidth)
        .attr('stroke-opacity', strokeOpacity)
        .attr('d', line);
    });

    // Add base calls and quality scores if available
    if (bases && peakLocations) {
      // Quality score bars
      const qvBarHeight = plotHeight * 0.15;
      // Move quality bars lower (closer to bottom)
      const qvBarY = yScale(yExtent[1] * 0.7); // Lowered from 0.8 to 0.7 for more space

      // Add tooltip div (one per plot, absolute position, hidden by default)
      let tooltip = d3.select(el).select('.chromatogram-tooltip');
      if (tooltip.empty()) {
        // Use the same dark tooltip styling as the sequence viewer for consistency
        tooltip = d3.select(el)
          .append('div')
          .attr('class', 'chromatogram-tooltip')
          .style('position', 'absolute')
          .style('pointer-events', 'none')
          .style('background', '#333')
          .style('color', 'white')
          .style('padding', '6px 10px')
          .style('border-radius', '4px')
          .style('font-size', '13px')
          .style('opacity', 0.95)
          .style('z-index', 1000)
          .style('display', 'none');
      }

      svg.selectAll('.qv-bar')
        .data(peakLocations)
        .enter()
        .append('line')
        .attr('class', 'qv-bar')
        .attr('x1', d => xScale(d))
        .attr('x2', d => xScale(d))
        .attr('y1', qvBarY)
        .attr('y2', (d, i) => {
          const qv = qualityScores[i];
          const height = qv ? (qv / 60 * qvBarHeight) : 0;
          return qvBarY - height;
        })
        .attr('stroke', (d, i) => {
          const qv = qualityScores[i];
          if (!qv) return QV_COLORS.fail;
          return qv >= 20 ? QV_COLORS.pass :
                 qv >= 15 ? QV_COLORS.suspect :
                           QV_COLORS.fail;
        })
        .attr('stroke-width', 4)
        .on('mouseover', function(event, d) {
          const i = peakLocations.indexOf(d);
          const qv = qualityScores[i];
          const base = bases && bases[i] ? bases[i] : '';
          tooltip.style('display', 'block')
            .html(`${base}:${i + 1} | QV:${qv}`);
          d3.select(this).attr('stroke-width', 6);
        })
        .on('mousemove', function(event) {
          // Position tooltip using mouse coordinates
          const x = event.pageX;
          const y = event.pageY;
          const tipNode = tooltip.node();
          const tipW = tipNode ? tipNode.offsetWidth : 0;
          const tipH = tipNode ? tipNode.offsetHeight : 0;
          
          // Position above the cursor
          tooltip
            .style('left', (x - tipW/2 - 250) + 'px')
            .style('top', (y - tipH - 0) + 'px');
        })
        .on('mouseout', function() {
          tooltip.style('display', 'none');
          d3.select(this).attr('stroke-width', 4);
        });

      // Base letters
      // Move bases higher (closer to top)
      const baseY = yScale(yExtent[1] * 0.95); // Already near top, but can be adjusted if needed
      // Bind bases together with their index so handlers can access quality & position
      const baseData = bases.map((b, i) => ({ base: b, idx: i }));
      svg.selectAll('.base-text')
        .data(baseData)
        .enter()
        .append('text')
        .attr('class', 'base-text')
        .attr('x', d => xScale(peakLocations[d.idx]))
        .attr('y', baseY)
        .attr('text-anchor', 'middle')
        .attr('font-family', 'monospace')
        .attr('font-size', '12px')
        .attr('font-weight', 'normal')
        .attr('fill', d => BASE_COLORS[d.base])
        .text(d => d.base)
        .style('cursor', 'default')
        .on('mouseover', function(event, d) {
          // Show tooltip with same style/content as sequence viewer
          const i = d.idx;
          const qv = qualityScores && qualityScores[i] ? qualityScores[i] : 'NA';
          const pos = i + 1;
          tooltip.style('display', 'block')
            .html(`${d.base}:${pos} | QV:${qv}`);
          d3.select(this).attr('font-weight', 'bold');
          //d3.select(this).attr('background-color', '#c3c3c3ff');
        })
        .on('mousemove', function(event) {
          // Position tooltip using mouse coordinates
          const x = event.pageX;
          const y = event.pageY;
          const tipNode = tooltip.node();
          const tipW = tipNode ? tipNode.offsetWidth : 0;
          const tipH = tipNode ? tipNode.offsetHeight : 0;
          
          // Position above the cursor
          tooltip
            .style('left', (x - tipW/2 - 250) + 'px')
            .style('top', (y - tipH - 0) + 'px');
        })
        .on('mouseout', function() {
          tooltip.style('display', 'none');
          d3.select(this).attr('font-weight', 'normal');
        });
    }

    // Add axes
    // Create tick values every 50 bases using peak locations
    const createTickValues = () => {
      if (!peakLocations || peakLocations.length === 0) return [];
      
      // Create mapping of base indices to x positions
      const baseToPos = new Map();
      peakLocations.forEach((pos, idx) => {
        baseToPos.set(idx + 1, pos);
      });
      
      // Create ticks every 50 bases
      const maxBases = peakLocations.length;
      const ticks = [];
      for (let i = 0; i <= maxBases; i += 50) {
        if (baseToPos.has(i)) {
          ticks.push(baseToPos.get(i));
        }
      }
      return ticks;
    };

    const xAxis = d3.axisBottom(xScale)
      .tickValues(createTickValues())
      .tickSizeOuter(0)
      .tickFormat((d) => {
        // Find the base index for this position
        const idx = peakLocations.indexOf(d);
        return idx >= 0 ? (idx + 1).toString() : '';
      });

    const yAxis = d3.axisLeft(yScale);

    svg.append('g')
      .attr('transform', `translate(0,${plotHeight})`)
      .attr('class', 'x-axis')
      .call(xAxis);

    svg.append('g')
      .call(yAxis);
  }

  // Register with Shiny if available
  if (window.Shiny) {
    const binding = new Shiny.OutputBinding();
    
    $.extend(binding, {
      find: function(scope) {
        return $(scope).find('.chromatogram-output');
      },
      renderValue: function(el, data) {
        renderChromatogram(el, data);
      },
      resize: function(el, width, height) {
        // Re-render on resize
        const data = $(el).data('chromatogram-data');
        if (data) renderChromatogram(el, data);
      }
    });

    Shiny.outputBindings.register(binding, 'shiny.chromatogramOutput');
  }
})();