package ca.mcgill.genome.bamtools.mutations;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.awt.Shape;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.jfree.chart.JFreeChart;
import org.jfree.chart.LegendItem;
import org.jfree.chart.LegendItemCollection;
import org.jfree.chart.LegendItemSource;
import org.jfree.chart.axis.AxisSpace;
import org.jfree.chart.axis.AxisState;
import org.jfree.chart.axis.CategoryAxis;
import org.jfree.chart.axis.CategoryLabelPosition;
import org.jfree.chart.axis.CategoryTick;
import org.jfree.chart.entity.CategoryLabelEntity;
import org.jfree.chart.entity.EntityCollection;
import org.jfree.chart.plot.Plot;
import org.jfree.chart.plot.PlotRenderingInfo;
import org.jfree.chart.title.LegendTitle;
import org.jfree.text.TextBlock;
import org.jfree.ui.RectangleAnchor;
import org.jfree.ui.RectangleEdge;

public class CovariateAxis extends CategoryAxis {
	private static final long serialVersionUID = 4219454789947501981L;
	private double markerSize = 15;
	private double extraMargin = markerSize/2;

	private final Map<String, Map<String, Integer>> covariates;
	
	public CovariateAxis(Map<String, Map<String, Integer>> covariates) {
		this.covariates = covariates;
	}

	public void addLegends(JFreeChart jFreeChart) {
		LegendTitle sex = new LegendTitle(new LegendItemSource() {
			@Override
			public LegendItemCollection getLegendItems() {
				LegendItemCollection lic = new LegendItemCollection();
			    Shape shape = new Rectangle(10, 10);
			    lic.add(new LegendItem("male", null, null, null, shape, Color.BLUE));
			    lic.add(new LegendItem("female", null, null, null, shape, Color.PINK));
				return lic;
			}
		});
		sex.setPosition(RectangleEdge.BOTTOM);

		LegendTitle age = new LegendTitle(new LegendItemSource() {
			@Override
			public LegendItemCollection getLegendItems() {
				LegendItemCollection lic = new LegendItemCollection();
			    Shape shape = new Rectangle(10, 10);
			    lic.add(new LegendItem("0-50", null, null, null, shape, Color.YELLOW));
			    lic.add(new LegendItem("51-60", null, null, null, shape, Color.ORANGE));
			    lic.add(new LegendItem("61-70", null, null, null, shape, Color.CYAN));
			    lic.add(new LegendItem(">70", null, null, null, shape, Color.RED));
				return lic;
			}
		});
		age.setPosition(RectangleEdge.BOTTOM);

		LegendTitle grade = new LegendTitle(new LegendItemSource() {
			@Override
			public LegendItemCollection getLegendItems() {
				LegendItemCollection lic = new LegendItemCollection();
			    Shape shape = new Rectangle(10, 10);
			    lic.add(new LegendItem("Grade 1", null, null, null, shape, new Color(0xADD8E6)));
			    lic.add(new LegendItem("Grade 2", null, null, null, shape, new Color(0x008000)));
			    lic.add(new LegendItem("Grade 3", null, null, null, shape, new Color(0xFFA500)));
			    lic.add(new LegendItem("Grade 4", null, null, null, shape, new Color(0x800000)));
				return lic;
			}
		});
		grade.setPosition(RectangleEdge.BOTTOM);
		
		//Order is important
		jFreeChart.addLegend(grade);
		jFreeChart.addLegend(age);
		jFreeChart.addLegend(sex);
	}
	@Override
	public AxisSpace reserveSpace(Graphics2D g2, Plot plot, Rectangle2D plotArea, RectangleEdge edge, AxisSpace space) {
		AxisSpace newSpace = super.reserveSpace(g2, plot, plotArea, edge, space);
		String anySample = covariates.keySet().iterator().next();
		if(RectangleEdge.isTopOrBottom(edge)) {
			newSpace.add(markerSize*covariates.get(anySample).size()+extraMargin, edge);
		}
		else {
			if((newSpace.getLeft()+newSpace.getRight()) < markerSize*covariates.get(anySample).size()) {
				newSpace.add(markerSize*covariates.get(anySample).size(), edge);
			}
		}
		return newSpace;
	}

	@Override
    @SuppressWarnings("rawtypes")
	protected AxisState drawCategoryLabels(Graphics2D g2, Rectangle2D plotArea,
			Rectangle2D dataArea, RectangleEdge edge, AxisState state,
			PlotRenderingInfo plotState) {
		AxisState newState = super.drawCategoryLabels(g2, plotArea, dataArea, edge, state, plotState);
		
		if (isTickLabelsVisible()) {
			List ticks =  state.getTicks();

            int categoryIndex = 0;
            Iterator iterator = ticks.iterator();
            while (iterator.hasNext()) {
                CategoryTick tick = (CategoryTick) iterator.next();
                String sampleName  = tick.getCategory().toString();

                CategoryLabelPosition position = getCategoryLabelPositions().getLabelPosition(edge);
                double x0 = 0.0;
                double x1 = 0.0;
                double y0 = 0.0;
                double y1 = 0.0;
                if (edge == RectangleEdge.TOP) {
                    x0 = getCategoryStart(categoryIndex, ticks.size(), dataArea, edge);
                    x1 = getCategoryEnd(categoryIndex, ticks.size(), dataArea, edge);
                    y1 = state.getCursor() - getCategoryLabelPositionOffset();
                    y0 = y1 - state.getMax();
                }
                else if (edge == RectangleEdge.BOTTOM) {
                    x0 = getCategoryStart(categoryIndex, ticks.size(), dataArea, edge);
                    x1 = getCategoryEnd(categoryIndex, ticks.size(), dataArea, edge);
                    y0 = state.getCursor() + getCategoryLabelPositionOffset();
                    y1 = y0 + state.getMax();
                }
                else if (edge == RectangleEdge.LEFT) {
                    y0 = getCategoryStart(categoryIndex, ticks.size(), dataArea, edge);
                    y1 = getCategoryEnd(categoryIndex, ticks.size(), dataArea, edge);
                    x1 = state.getCursor() - getCategoryLabelPositionOffset();
                    x0 = x1 - state.getMax();
                }
                else if (edge == RectangleEdge.RIGHT) {
                    y0 = getCategoryStart(categoryIndex, ticks.size(), dataArea, edge);
                    y1 = getCategoryEnd(categoryIndex, ticks.size(), dataArea, edge);
                    x0 = state.getCursor() + getCategoryLabelPositionOffset();
                    x1 = x0 - state.getMax();
                }
                Rectangle2D area = new Rectangle2D.Double(x0, y0, (x1 - x0), (y1 - y0));
                Point2D anchorPoint = RectangleAnchor.coordinates(area, position.getCategoryAnchor());
                TextBlock block = tick.getLabel();
                Shape bounds = block.calculateBounds(g2,(float) anchorPoint.getX(), (float) anchorPoint.getY(),position.getLabelAnchor(), (float) anchorPoint.getX(),(float) anchorPoint.getY(), position.getAngle());
                
        		Color c = g2.getColor();
//        		g2.setColor(Color.BLACK);
//                g2.drawRect(bounds.getBounds().x, bounds.getBounds().y,bounds.getBounds().width, bounds.getBounds().height);
                int xCenter = bounds.getBounds().x + bounds.getBounds().width/2;
                int offset=0;
                if(sampleName != null && sampleName.length() > 0) {
	        		g2.setColor(Color.WHITE);
	                if(covariates.get(sampleName).get("sex").equals(1)) {
		        		g2.setColor(Color.BLUE);
	                }
	                else if(covariates.get(sampleName).get("sex").equals(2)) {
		        		g2.setColor(Color.PINK);
	                }
	                g2.fillRect((int)(xCenter-(markerSize/2)), (int)(bounds.getBounds().y+offset), (int)markerSize, (int)markerSize);
	        		g2.setColor(Color.WHITE);
	                g2.drawRect((int)(xCenter-(markerSize/2)), (int)(bounds.getBounds().y+offset), (int)markerSize, (int)markerSize);
	                offset += markerSize;
	                
	                int age = covariates.get(sampleName).get("age").intValue();
	                if(age == 0) {
	                	g2.setColor(Color.WHITE);
	                }
	                else if(age < 51) {
		        		g2.setColor(Color.YELLOW);
	                }
	                else if(age < 61) {
		        		g2.setColor(Color.ORANGE);
	                }
	                else if(age < 71) {
		        		g2.setColor(Color.CYAN);
	                }
	                else {
		        		g2.setColor(Color.RED);
	                }
	                
	                g2.fillRect((int)(xCenter-(markerSize/2)), (int)(bounds.getBounds().y+offset), (int)markerSize, (int)markerSize);
	        		g2.setColor(Color.WHITE);
	                g2.drawRect((int)(xCenter-(markerSize/2)), (int)(bounds.getBounds().y+offset), (int)markerSize, (int)markerSize);
	                offset += markerSize;
	                int grade = covariates.get(sampleName).get("grade").intValue();
                	g2.setColor(Color.WHITE);
	                if(grade == 1) {
		        		g2.setColor(new Color(0xADD8E6));
	                }
	                else if(grade == 2) {
		        		g2.setColor(new Color(0x008000));
	                }
	                else if(grade == 3) {
		        		g2.setColor(new Color(0xFFA500));
	                }
	                else if(grade == 4) {
		        		g2.setColor(new Color(0x800000));
	                }
	                g2.fillRect((int)(xCenter-(markerSize/2)), (int)(bounds.getBounds().y+offset), (int)markerSize, (int)markerSize);
	        		g2.setColor(Color.WHITE);
	                g2.drawRect((int)(xCenter-(markerSize/2)), (int)(bounds.getBounds().y+offset), (int)markerSize, (int)markerSize);
	                offset += markerSize;
/*	                
	                for(String covName : covariates.get(sampleName).keySet()) {
	                	
		        		g2.setColor(Color.RED);
		                g2.fillRect((int)(bounds.getBounds().x+(markerSize/2)), (int)(bounds.getBounds().y+offset), (int)markerSize, (int)markerSize);
		        		g2.setColor(Color.WHITE);
		                g2.drawRect((int)(bounds.getBounds().x+(markerSize/2)), (int)(bounds.getBounds().y+offset), (int)markerSize, (int)markerSize);
		                offset += markerSize;
	                }
*/
                }
        		g2.setColor(c);
                
                
                bounds = new Rectangle2D.Double(bounds.getBounds().x, bounds.getBounds().y, bounds.getBounds().width, bounds.getBounds().height+20);
                if (plotState != null && plotState.getOwner() != null) {
                    EntityCollection entities
                            = plotState.getOwner().getEntityCollection();
                    if (entities != null) {
                        String tooltip = getCategoryLabelToolTip(
                                tick.getCategory());
                        entities.add(new CategoryLabelEntity(tick.getCategory(),
                                bounds, tooltip, null));
                    }
                }
                categoryIndex++;
            }
		}
        return newState;
	}
	
	
}
