// Sequence viewer component for browser-side sequence rendering
(function() {
  const QV_COLORS = {
    high: '#4CAF50',    // Green (QV >= 20)
    medium: '#FFC107',  // Amber (QV >= 10)
    low: '#F44336'      // Red (QV < 10)
  };

  function renderSequence(el, data) {
    console.log('renderSequence called with data:', data);
    if (!data || !data.bases || !data.qscores) {
      console.warn('Invalid sequence data:', data);
      return;
    }
    const { bases, qscores, crlStart, crlEnd } = data;
    if (bases.length !== qscores.length) {
      console.warn('Bases and qscores length mismatch:', bases.length, qscores.length);
      return;
    }

    // Clear existing content
    el.innerHTML = '';
    
    // Create container
    const container = document.createElement('div');
    container.style.cssText = `
      white-space: pre-wrap;
      word-break: break-all;
      min-height: 100px;
      max-height: 600px;
      overflow-y: auto;
      border: 1px solid #ddd;
      padding: 14px;
      font-size: 11px;
      font-family: monospace;
      line-height: 1.5;
    `;

    // Create base spans with quality coloring
    for (let i = 0; i < bases.length; i++) {
      const span = document.createElement('span');
      const qv = qscores[i];
      const isInCRL = i >= crlStart && i <= crlEnd;
      
      // Base styling
      span.textContent = bases[i];
      span.style.color = qv >= 20 ? QV_COLORS.high :
                        qv >= 10 ? QV_COLORS.medium :
                                 QV_COLORS.low;
      span.style.backgroundColor = !isInCRL ? '#c3c3c3ff' : 'transparent';
      span.style.padding = '1px 0';
      span.style.cursor = 'default';
      
      // Tooltip data
      span.dataset.position = i + 1;
      span.dataset.qv = Math.round(qv);
      
      // Add hover effect and tooltip
      span.onmouseover = function(e) {
        this.style.backgroundColor = '#dac586';
        showTooltip(e, `${bases[i]}:${i + 1} | QV:${Math.round(qv)}`);
      };
      span.onmouseout = function() {
        this.style.backgroundColor = !isInCRL ? '#c3c3c3ff' : 'transparent';
        hideTooltip();
      };
      span.onmousemove = updateTooltipPosition;
      
      container.appendChild(span);
      
      // Add spacing every 10 bases and line break every 100 bases
      if (i < bases.length - 1) {
        if ((i + 1) % 100 === 0) {
          container.appendChild(document.createElement('br'));
        } else if ((i + 1) % 10 === 0) {
          container.appendChild(document.createTextNode('\u00A0')); // &nbsp;
        }
      }
    }
    
    el.appendChild(container);
  }

  // Tooltip handling
  let tooltip = null;

  function createTooltip() {
    tooltip = document.createElement('div');
    tooltip.style.cssText = `
      position: fixed;
      background: #333;
      color: white;
      padding: 4px 8px;
      border-radius: 4px;
      font-size: 13px;
      pointer-events: none;
      opacity: 0.9;
      z-index: 1000;
      display: none;
    `;
    document.body.appendChild(tooltip);
  }

  function showTooltip(event, text) {
    if (!tooltip) createTooltip();
    tooltip.textContent = text;
    tooltip.style.display = 'block';
    updateTooltipPosition(event);
  }
  
  // Initialize Shiny binding
  const sequenceBinding = new Shiny.OutputBinding();

  $.extend(sequenceBinding, {
    find: function(scope) {
      return $(scope).find('.sequence-viewer-container');
    },

    renderValue: function(el, data) {
      if (!data) {
        console.warn('No sequence data provided');
        return;
      }
      console.log('Rendering sequence data:', data);
      renderSequence(el, data);
    }
  });

  // Register the binding
  Shiny.outputBindings.register(sequenceBinding, 'shiny.sequenceViewer');

  function updateTooltipPosition(event) {
    if (!tooltip) return;
    const x = event.pageX + 15;
    const y = event.pageY - 15;
    tooltip.style.left = x + 'px';
    tooltip.style.top = y + 'px';
  }

  function hideTooltip() {
    if (tooltip) tooltip.style.display = 'none';
  }

  // Handle window resize to adjust container size
  window.addEventListener('resize', function() {
    const containers = document.querySelectorAll('.sequence-output > div');
    containers.forEach(container => {
      // Adjust max-height based on window height
      const maxHeight = Math.max(100, Math.min(600, window.innerHeight * 0.6));
      container.style.maxHeight = maxHeight + 'px';
    });
  });
})();