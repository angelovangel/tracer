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
    
    // Create container with optimized styles
    const container = document.createElement('div');
    Object.assign(container.style, {
      whiteSpace: 'pre-wrap',
      wordBreak: 'break-all',
      minHeight: '100px',
      maxHeight: '600px',
      overflowY: 'auto',
      border: '1px solid #ddd',
      padding: '14px',
      fontSize: '13px',
      fontFamily: 'monospace',
      lineHeight: '1.5'
    });

    // Add event delegation for better performance
    container.addEventListener('mouseover', function(e) {
      const span = e.target;
      if (span.tagName === 'SPAN') {
        span.style.backgroundColor = '#dac586';
        showTooltip(e, `${span.textContent}:${+span.dataset.i + 1} | QV:${span.dataset.qv}`);
      }
    });

    container.addEventListener('mouseout', function(e) {
      const span = e.target;
      if (span.tagName === 'SPAN') {
        span.style.backgroundColor = span.dataset.inCrl === 'false' ? '#c3c3c3ff' : 'transparent';
        hideTooltip();
      }
    });

    container.addEventListener('mousemove', function(e) {
      if (e.target.tagName === 'SPAN') {
        updateTooltipPosition(e);
      }
    });

    // Create all base spans at once using DocumentFragment
    const fragment = document.createDocumentFragment();
    const baseElements = new Array(bases.length);
    
    // Pre-create common elements for spacing
    const spaceNode = document.createTextNode('\u00A0');
    const brNode = document.createElement('br');
    
    // Pre-compute base styles for better performance
    const baseStyle = {
      padding: '1px 0',
      cursor: 'default'
    };
    
    // Use a single template for better performance
    const spanTemplate = document.createElement('span');
    Object.assign(spanTemplate.style, baseStyle);
    
    // Build all elements in memory first
    for (let i = 0; i < bases.length; i++) {
      const span = spanTemplate.cloneNode(false);
      const qv = qscores[i];
      const isInCRL = i >= crlStart && i <= crlEnd;
      
      // Set content and styling
      span.textContent = bases[i];
      span.style.color = qv >= 20 ? QV_COLORS.high :
                        qv >= 10 ? QV_COLORS.medium :
                                 QV_COLORS.low;
      span.style.backgroundColor = !isInCRL ? '#c3c3c3ff' : 'transparent';
      
      // Store minimal data needed for tooltips
      span.dataset.i = i;
      span.dataset.qv = qv;
      span.dataset.inCrl = isInCRL;
      
      // Add to fragment
      fragment.appendChild(span);
      baseElements[i] = span;
      
      // Add spacing using cached nodes
      if (i < bases.length - 1) {
        if ((i + 1) % 100 === 0) {
          fragment.appendChild(brNode.cloneNode());
        } else if ((i + 1) % 10 === 0) {
          fragment.appendChild(spaceNode.cloneNode());
        }
      }
    }
    
    // Append the fragment to container first
    container.appendChild(fragment);
    
    // Then append container to the element
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