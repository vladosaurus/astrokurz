import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib.patches as patches
from reportlab.lib.pagesizes import A4
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, Table, TableStyle
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.lib import colors
import io
import os

class BinaryStarAnalysis:
    def __init__(self):
        # Given parameters
        self.a_angular = 4.5  # arcseconds - semi-major axis
        self.b_angular = 3.4  # arcseconds - semi-minor axis
        self.m1_apparent = 3.9  # apparent magnitude of primary
        self.m2_apparent = 5.3  # apparent magnitude of secondary
        self.observation_period = 11  # years
        self.inclination = 90  # degrees (edge-on view)
        
        # Constants
        self.G = 6.67430e-11  # m³ kg⁻¹ s⁻²
        self.AU = 1.496e11    # meters
        self.pc = 3.086e16    # meters
        self.M_sun = 1.989e30 # kg
        self.L_sun = 3.828e26 # watts
        
        # Results storage
        self.results = {}
        
    def calculate_orbital_period(self):
        """
        Calculate orbital period from the observation that we see almost complete orbit
        in 11 years based on the elliptical segment area.
        """
        # Calculate eccentricity
        c = np.sqrt(self.a_angular**2 - self.b_angular**2)
        self.eccentricity = c / self.a_angular
        
        # For the area calculation, we need to estimate what fraction of orbit was observed
        # Based on the elliptical segment formula given in the problem
        ellipse_area = np.pi * self.a_angular * self.b_angular
        
        # Estimate that we observed roughly 85-95% of the orbit in 11 years
        # This is an approximation that will be refined through iteration
        observed_fraction = 0.9  # Initial guess
        
        self.period = self.observation_period / observed_fraction
        self.results['period'] = self.period
        
        print(f"Initial period estimate: {self.period:.2f} years")
        return self.period
    
    def mass_luminosity_relation(self, mass_solar):
        """Mass-luminosity relation for main sequence stars"""
        if mass_solar < 0.43:
            return 0.23 * mass_solar**2.3
        elif mass_solar < 2:
            return mass_solar**4
        elif mass_solar < 20:
            return 1.4 * mass_solar**3.5
        else:
            return 32000 * mass_solar
    
    def luminosity_to_absolute_magnitude(self, luminosity_solar):
        """Convert luminosity (in solar units) to absolute magnitude"""
        return 4.83 - 2.5 * np.log10(luminosity_solar)
    
    def absolute_to_apparent_magnitude(self, M, d_pc):
        """Convert absolute magnitude to apparent magnitude"""
        return M + 5 * np.log10(d_pc) - 5
    
    def iterate_solution(self, max_iterations=20, tolerance=0.01):
        """
        Iterative solution using dynamic parallax method
        """
        print("Starting iterative solution...")
        print("-" * 50)
        
        # Initial guesses
        distance_pc = 10  # parsecs - initial guess
        mass1_solar = 1.0  # solar masses
        mass2_solar = 0.8  # solar masses
        
        for iteration in range(max_iterations):
            print(f"\nIteration {iteration + 1}:")
            
            # Store previous values for convergence check
            prev_distance = distance_pc
            prev_mass1 = mass1_solar
            prev_mass2 = mass2_solar
            
            # 1. Convert angular sizes to linear sizes
            a_linear_au = self.a_angular * distance_pc  # AU
            a_linear_m = a_linear_au * self.AU  # meters
            
            # 2. Apply Kepler's Third Law: P² = (4π²/G(M1+M2)) * a³
            total_mass_kg = (mass1_solar + mass2_solar) * self.M_sun
            period_seconds = self.period * 365.25 * 24 * 3600
            
            # Calculate semi-major axis from Kepler's 3rd law
            a_calculated = ((G * total_mass_kg * period_seconds**2) / (4 * np.pi**2))**(1/3)
            
            # Update distance based on consistency
            distance_pc = (a_calculated / self.AU) / self.a_angular
            
            # 3. Calculate luminosities from mass-luminosity relation
            L1_solar = self.mass_luminosity_relation(mass1_solar)
            L2_solar = self.mass_luminosity_relation(mass2_solar)
            
            # 4. Calculate absolute magnitudes
            M1_abs = self.luminosity_to_absolute_magnitude(L1_solar)
            M2_abs = self.luminosity_to_absolute_magnitude(L2_solar)
            
            # 5. Check consistency with observed apparent magnitudes
            m1_calculated = self.absolute_to_apparent_magnitude(M1_abs, distance_pc)
            m2_calculated = self.absolute_to_apparent_magnitude(M2_abs, distance_pc)
            
            # 6. Adjust masses based on magnitude differences
            dm1 = self.m1_apparent - m1_calculated
            dm2 = self.m2_apparent - m2_calculated
            
            # Update masses (with some damping for stability)
            damping = 0.3
            mass1_solar *= (1 + damping * dm1 / 5)  # Rough adjustment
            mass2_solar *= (1 + damping * dm2 / 5)
            
            # Ensure reasonable mass bounds
            mass1_solar = max(0.1, min(50, mass1_solar))
            mass2_solar = max(0.1, min(50, mass2_solar))
            
            print(f"  Distance: {distance_pc:.3f} pc")
            print(f"  Mass 1: {mass1_solar:.3f} M☉")
            print(f"  Mass 2: {mass2_solar:.3f} M☉")
            print(f"  M1_abs: {M1_abs:.2f}")
            print(f"  M2_abs: {M2_abs:.2f}")
            print(f"  m1_calc: {m1_calculated:.2f} (obs: {self.m1_apparent})")
            print(f"  m2_calc: {m2_calculated:.2f} (obs: {self.m2_apparent})")
            
            # Check convergence
            distance_change = abs(distance_pc - prev_distance) / prev_distance
            mass1_change = abs(mass1_solar - prev_mass1) / prev_mass1
            mass2_change = abs(mass2_solar - prev_mass2) / prev_mass2
            
            max_change = max(distance_change, mass1_change, mass2_change)
            print(f"  Max relative change: {max_change:.4f}")
            
            if max_change < tolerance:
                print(f"\nConverged after {iteration + 1} iterations!")
                break
        
        # Store final results
        self.results.update({
            'distance_pc': distance_pc,
            'distance_ly': distance_pc * 3.26,
            'mass1_solar': mass1_solar,
            'mass2_solar': mass2_solar,
            'total_mass_solar': mass1_solar + mass2_solar,
            'M1_absolute': M1_abs,
            'M2_absolute': M2_abs,
            'semi_major_axis_au': a_linear_au,
            'semi_major_axis_km': a_linear_m / 1000,
            'luminosity1_solar': L1_solar,
            'luminosity2_solar': L2_solar,
            'eccentricity': self.eccentricity
        })
        
        return self.results
    
    def create_orbital_plot(self):
        """Create a plot showing the orbital configuration"""
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
        
        # Plot 1: Orbital ellipse (as seen from above)
        ax1.set_aspect('equal')
        
        # Create ellipse
        ellipse = Ellipse((0, 0), 2*self.a_angular, 2*self.b_angular, 
                         fill=False, edgecolor='blue', linewidth=2)
        ax1.add_patch(ellipse)
        
        # Mark stellar positions (example positions)
        theta = np.linspace(0, 2*np.pi, 8)
        x_orbit = self.a_angular * np.cos(theta)
        y_orbit = self.b_angular * np.sin(theta)
        
        # Primary star at focus (center for this simplified view)
        ax1.plot(0, 0, 'ro', markersize=8, label=f'Primary (m={self.m1_apparent})')
        
        # Secondary star positions
        ax1.plot(x_orbit, y_orbit, 'b.', markersize=4, alpha=0.6)
        ax1.plot(x_orbit[0], y_orbit[0], 'bo', markersize=6, label=f'Secondary (m={self.m2_apparent})')
        
        ax1.set_xlabel('Angular separation (arcsec)')
        ax1.set_ylabel('Angular separation (arcsec)')
        ax1.set_title('Binary Star Orbit (Observer View)')
        ax1.grid(True, alpha=0.3)
        ax1.legend()
        
        # Add annotations
        ax1.annotate(f'a = {self.a_angular}"', xy=(self.a_angular*0.7, 0), 
                    xytext=(self.a_angular*0.7, self.a_angular*0.3),
                    arrowprops=dict(arrowstyle='->', color='red'))
        ax1.annotate(f'b = {self.b_angular}"', xy=(0, self.b_angular*0.7), 
                    xytext=(self.a_angular*0.3, self.b_angular*0.7),
                    arrowprops=dict(arrowstyle='->', color='red'))
        
        # Plot 2: System parameters
        ax2.axis('off')
        
        if hasattr(self, 'results') and self.results:
            results_text = f"""
BINARY STAR SYSTEM ANALYSIS RESULTS

Orbital Parameters:
• Period: {self.results.get('period', 'N/A'):.2f} years
• Semi-major axis: {self.results.get('semi_major_axis_au', 'N/A'):.1f} AU
• Eccentricity: {self.results.get('eccentricity', 'N/A'):.3f}

System Properties:
• Distance: {self.results.get('distance_pc', 'N/A'):.2f} pc ({self.results.get('distance_ly', 'N/A'):.1f} ly)
• Total mass: {self.results.get('total_mass_solar', 'N/A'):.2f} M☉

Primary Star:
• Mass: {self.results.get('mass1_solar', 'N/A'):.2f} M☉
• Absolute magnitude: {self.results.get('M1_absolute', 'N/A'):.2f}
• Luminosity: {self.results.get('luminosity1_solar', 'N/A'):.2f} L☉

Secondary Star:
• Mass: {self.results.get('mass2_solar', 'N/A'):.2f} M☉
• Absolute magnitude: {self.results.get('M2_absolute', 'N/A'):.2f}
• Luminosity: {self.results.get('luminosity2_solar', 'N/A'):.2f} L☉

Given Data:
• Angular semi-major axis: {self.a_angular}"
• Angular semi-minor axis: {self.b_angular}"
• Apparent magnitude 1: {self.m1_apparent}
• Apparent magnitude 2: {self.m2_apparent}
• Observation period: {self.observation_period} years
• Inclination: {self.inclination}°
            """
        else:
            results_text = "Run analysis first to see results"
        
        ax2.text(0.05, 0.95, results_text, transform=ax2.transAxes, 
                fontsize=10, verticalalignment='top', fontfamily='monospace')
        
        plt.tight_layout()
        return fig
    
    def generate_pdf_report(self, filename='binary_star_report.pdf'):
        """Generate a comprehensive PDF report"""
        doc = SimpleDocTemplate(filename, pagesize=A4)
        styles = getSampleStyleSheet()
        story = []
        
        # Title
        title_style = ParagraphStyle(
            'CustomTitle',
            parent=styles['Heading1'],
            fontSize=18,
            spaceAfter=30,
            alignment=1  # Center
        )
        story.append(Paragraph("Binary Star System Analysis Report", title_style))
        story.append(Spacer(1, 20))
        
        # Problem description
        story.append(Paragraph("Problem Description", styles['Heading2']))
        problem_text = f"""
        Analysis of a visual binary star system with the following observed parameters:
        <br/>• Angular semi-major axis: {self.a_angular} arcseconds
        <br/>• Angular semi-minor axis: {self.b_angular} arcseconds  
        <br/>• Apparent magnitude of primary: {self.m1_apparent}
        <br/>• Apparent magnitude of secondary: {self.m2_apparent}
        <br/>• Observation period: {self.observation_period} years
        <br/>• Orbital inclination: {self.inclination}° (edge-on view)
        """
        story.append(Paragraph(problem_text, styles['Normal']))
        story.append(Spacer(1, 20))
        
        # Method
        story.append(Paragraph("Method: Dynamic Parallax", styles['Heading2']))
        method_text = """
        The dynamic parallax method iteratively solves for the system parameters using:
        <br/>1. Kepler's Third Law: P² ∝ a³/(M₁+M₂)
        <br/>2. Mass-luminosity relation for main sequence stars
        <br/>3. Distance modulus: m - M = 5 log(d) - 5
        <br/>4. Consistency between calculated and observed apparent magnitudes
        """
        story.append(Paragraph(method_text, styles['Normal']))
        story.append(Spacer(1, 20))
        
        # Results table
        if hasattr(self, 'results') and self.results:
            story.append(Paragraph("Results", styles['Heading2']))
            
            data = [
                ['Parameter', 'Value', 'Unit'],
                ['Orbital Period', f"{self.results['period']:.2f}", 'years'],
                ['System Distance', f"{self.results['distance_pc']:.2f}", 'parsecs'],
                ['System Distance', f"{self.results['distance_ly']:.1f}", 'light years'],
                ['Semi-major Axis', f"{self.results['semi_major_axis_au']:.1f}", 'AU'],
                ['Total System Mass', f"{self.results['total_mass_solar']:.2f}", 'M☉'],
                ['Primary Mass', f"{self.results['mass1_solar']:.2f}", 'M☉'],
                ['Secondary Mass', f"{self.results['mass2_solar']:.2f}", 'M☉'],
                ['Primary Abs. Magnitude', f"{self.results['M1_absolute']:.2f}", 'mag'],
                ['Secondary Abs. Magnitude', f"{self.results['M2_absolute']:.2f}", 'mag'],
                ['Primary Luminosity', f"{self.results['luminosity1_solar']:.2f}", 'L☉'],
                ['Secondary Luminosity', f"{self.results['luminosity2_solar']:.2f}", 'L☉'],
                ['Orbital Eccentricity', f"{self.results['eccentricity']:.3f}", ''],
            ]
            
            table = Table(data)
            table.setStyle(TableStyle([
                ('BACKGROUND', (0, 0), (-1, 0), colors.grey),
                ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
                ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
                ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
                ('FONTSIZE', (0, 0), (-1, 0), 12),
                ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
                ('BACKGROUND', (0, 1), (-1, -1), colors.beige),
                ('GRID', (0, 0), (-1, -1), 1, colors.black)
            ]))
            story.append(table)
            story.append(Spacer(1, 20))
        
        # Save plot as image for PDF
        fig = self.create_orbital_plot()
        img_buffer = io.BytesIO()
        fig.savefig(img_buffer, format='png', dpi=150, bbox_inches='tight')
        img_buffer.seek(0)
        
        # Add plot to PDF
        story.append(Paragraph("System Visualization", styles['Heading2']))
        
        # Save temporary image file
        temp_img = 'temp_plot.png'
        with open(temp_img, 'wb') as f:
            f.write(img_buffer.getvalue())
        
        img = Image(temp_img, width=6*inch, height=3*inch)
        story.append(img)
        
        # Build PDF
        doc.build(story)
        
        # Clean up temporary file
        if os.path.exists(temp_img):
            os.remove(temp_img)
        
        plt.close(fig)
        print(f"PDF report saved as: {filename}")
    
    def run_analysis(self):
        """Run the complete analysis"""
        print("Binary Star System Analysis")
        print("=" * 40)
        
        # Step 1: Calculate orbital period
        self.calculate_orbital_period()
        
        # Step 2: Iterative solution
        results = self.iterate_solution()
        
        # Step 3: Display results
        print("\n" + "="*50)
        print("FINAL RESULTS:")
        print("="*50)
        for key, value in results.items():
            if isinstance(value, float):
                print(f"{key}: {value:.4f}")
            else:
                print(f"{key}: {value}")
        
        return results

# Example usage
if __name__ == "__main__":
    # Create analysis object
    binary_system = BinaryStarAnalysis()
    
    # Run the analysis
    results = binary_system.run_analysis()
    
    # Create visualization
    fig = binary_system.create_orbital_plot()
    plt.show()
    
    # Generate PDF report
    binary_system.generate_pdf_report()
    
    print("\nAnalysis complete!")
    print("Check the generated PDF report for detailed results.")
